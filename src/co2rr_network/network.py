"""Construct reaction networks from one formula-to-energy dictionary."""
import hashlib
import json
from .chemistry import che_shift, state_composition
from .io import finite
from .rules import enumerate_states, load_rules
from .schema import Reaction, ReactionNetwork, State, state_id


def _parts(formula):
    parts = tuple(p.strip() for p in formula.split("+"))
    state_id(parts)
    return parts


def _entry(value, scalar_kind):
    if not isinstance(value, dict):
        return {scalar_kind: finite(value, scalar_kind)}
    if not value or set(value) - {"free_energy", "delta_g", "barrier"}:
        raise ValueError("Energy values accept free_energy, delta_g and barrier only")
    return {key: finite(v, key) for key, v in value.items()}


def build_network(energies, rules=None, profile="direct", energy_type="free_energy",
                  missing_energy="skip", potential=0.0, ph=0.0, temperature=298.15,
                  potential_scale="RHE", include_reverse=False, aliases=None, metadata=None):
    """Build from dict[str, float | dict[str, float]].

    Species/state keys carry formation G. Reaction keys, e.g. 'CO+CO -> COCO',
    carry delta_g and/or barrier; a scalar reaction value defaults to delta_g.
    Set energy_type='barrier' for a scalar barrier dictionary. Target-only edge
    values are accepted only if the supplied rules have one unique predecessor.
    Missing absolute state energies remain None rather than being invented.
    """
    if not isinstance(energies, dict) or not energies:
        raise ValueError("energies must be a nonempty formula-to-energy dict")
    if energy_type not in ("free_energy", "delta_g", "barrier"):
        raise ValueError("energy_type must be free_energy, delta_g or barrier")
    if missing_energy not in ("error", "skip"):
        raise ValueError("missing_energy must be error or skip")
    che_shift({}, potential, ph, temperature, potential_scale)
    state_values, targets, edges, normalized = {}, {}, {}, {}
    for formula, value in energies.items():
        if not isinstance(formula, str) or not formula.strip():
            raise ValueError("Energy dictionary keys must be nonempty formula strings")
        formula = formula.replace("→", "->")
        if "->" in formula:
            if formula.count("->") != 1:
                raise ValueError("A reaction key needs one arrow: " + formula)
            left, right = (_parts(s.strip()) for s in formula.split("->"))
            key = (state_id(left), state_id(right))
            if key in edges:
                raise ValueError("Duplicate normalized reaction key: " + formula)
            entry = _entry(value, "delta_g" if energy_type == "free_energy" else energy_type)
            if "free_energy" in entry:
                raise ValueError("A reaction key uses delta_g and/or barrier, not free_energy")
            edges[key] = entry
            normalized[key[0] + " -> " + key[1]] = entry
        else:
            key = state_id(_parts(formula))
            entry = _entry(value, energy_type)
            normalized[key] = entry
            if key in state_values or key in targets:
                raise ValueError("Duplicate normalized formula: " + key)
            if "free_energy" in entry:
                state_values[key] = entry["free_energy"]
            edge_values = {k: v for k, v in entry.items() if k != "free_energy"}
            if edge_values:
                targets[key] = edge_values
    if rules is None and edges and not state_values and not targets:
        inferred = [{"source": s.split("+"), "target": t.split("+"), "kind": "reaction"}
                    for s, t in edges]
        catalog = load_rules(data=inferred)
    else:
        catalog = load_rules(profile, data=rules) if rules is not None else load_rules(profile)
        known = {(state_id(r["source"]), state_id(r["target"])) for r in catalog["reactions"]}
        additional = [{"source": s.split("+"), "target": t.split("+"), "kind": "reaction"}
                      for s, t in edges if (s, t) not in known]
        if additional:
            catalog = load_rules(data=dict(catalog, reactions=catalog["reactions"] + additional))
    for target, values in targets.items():
        sources = {state_id(r["source"]) for r in catalog["reactions"] if state_id(r["target"]) == target}
        if len(sources) != 1:
            raise ValueError("{} has {} candidate predecessors; use a 'source -> target' reaction key"
                             .format(target, len(sources)))
        key = (next(iter(sources)), target)
        if key in edges:
            raise ValueError("Energy supplied twice for {} -> {}".format(*key))
        edges[key] = values
    formulas, aliases = catalog.get("formulas", {}), aliases or {}
    states, missing = {}, []
    for parts in enumerate_states(catalog):
        sid = state_id(parts)
        lookup_parts = tuple(p if p in state_values else aliases.get(p, p) for p in parts)
        lookup = sid if sid in state_values else state_id(lookup_parts)
        if lookup in state_values:
            g, model = state_values[lookup], "explicit_state"
        elif all(p in state_values for p in lookup_parts):
            g, model = sum(state_values[p] for p in lookup_parts), "independent_adsorbates"
        else:
            g, model = None, "not_supplied"
            missing.append(sid)
        counts = state_composition(parts, formulas)
        shift = che_shift(counts, potential, ph, temperature, potential_scale)
        states[sid] = State(sid, parts, None if g is None else g + shift,
                            any(p.endswith("_d") for p in parts),
                            {"composition": counts, "energy_model": model, "che_shift": shift})
    reactions = []
    skipped = []
    for rule in catalog["reactions"]:
        source, target = state_id(rule["source"]), state_id(rule["target"])
        known = edges.get((source, target), {})
        left, right = states[source], states[target]
        derived = (right.free_energy - left.free_energy
                   if right.free_energy is not None and left.free_energy is not None else None)
        dg = known.get("delta_g")
        if dg is not None:
            dg += right.metadata["che_shift"] - left.metadata["che_shift"]
            if derived is not None and abs(derived - dg) > 1e-7:
                raise ValueError("State energies and supplied delta_g disagree: {} -> {}".format(source, target))
        else:
            dg = derived
        barrier = known.get("barrier")
        if dg is None and barrier is None:
            skipped.append(source + " -> " + target)
            continue
        if barrier is not None and barrier + 1e-7 < max(0.0, dg if dg is not None else 0.0):
            raise ValueError("Barrier below a reaction endpoint: {} -> {}".format(source, target))
        reactions.append(Reaction("r{:05d}".format(len(reactions)), source, target, rule["kind"],
                                  dg, barrier, dict(rule.get("metadata", {}), **rule["balance"])))
    if skipped and missing_energy == "error":
        raise ValueError("Missing energy for {} reactions, e.g. {}".format(len(skipped), skipped[0]))
    if not reactions:
        raise ValueError("No reactions have usable energy data")
    if include_reverse:
        existing = {(r.source, r.target, r.kind) for r in reactions}
        for r in list(reactions):
            if (r.target, r.source, r.kind) in existing:
                continue
            if r.delta_g is None:
                raise ValueError("Reverse reactions require delta_g: " + r.id)
            barrier = None if r.barrier is None else max(0.0, r.barrier - r.delta_g)
            reactions.append(Reaction("r{:05d}".format(len(reactions)), r.target, r.source, r.kind,
                                      -r.delta_g, barrier, {"reverse_of": r.id}))
            existing.add((r.target, r.source, r.kind))
    used = {sid for r in reactions for sid in (r.source, r.target)}
    digest = hashlib.sha256(json.dumps(normalized, sort_keys=True, allow_nan=False).encode("utf-8")).hexdigest()
    return ReactionNetwork([s for s in states.values() if s.id in used], reactions,
                           dict(metadata or {}, conditions={"potential": potential, "ph": ph,
                               "temperature": temperature, "potential_scale": potential_scale},
                               energy_data_sha256=digest, rules_sha256=catalog["catalog_sha256"],
                               missing_state_energies=missing, skipped_reactions=skipped))
