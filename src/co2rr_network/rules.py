"""Curated CO2RR catalog plus a user-defined rule schema."""
import json
import hashlib
from pathlib import Path
from .chemistry import balanced_step
from .io import read_config, sha256
from .schema import state_id


def load_rules(profile="direct", path=None, data=None):
    if data is not None:
        data = dict(data) if isinstance(data, dict) else {"reactions": data}
        catalog_path = None
    elif path:
        data = read_config(path)
        catalog_path = Path(path)
    else:
        if profile not in ("direct", "out"):
            raise ValueError("Unknown rule profile: {}".format(profile))
        catalog_path = Path(__file__).parent / "resources" / (profile + "_rules.json")
        with catalog_path.open(encoding="utf-8") as stream:
            data = json.load(stream)
    formulas = data.get("formulas", {})
    seen, rules = set(), []
    for rule in data["reactions"]:
        if not isinstance(rule["source"], list) or not isinstance(rule["target"], list):
            raise ValueError("Rule source and target must be arrays of species names")
        source, target = tuple(rule["source"]), tuple(rule["target"])
        key = (source, target, rule["kind"])
        state_id(source)
        state_id(target)
        balance = balanced_step(source, target, formulas)
        if key not in seen and source != target:
            rules.append(dict(rule, source=list(source), target=list(target), balance=balance))
            seen.add(key)
    digest = sha256(catalog_path) if catalog_path else hashlib.sha256(
        json.dumps(data, sort_keys=True).encode("utf-8")).hexdigest()
    return dict(data, reactions=rules, catalog_sha256=digest)


def enumerate_states(rules):
    return sorted({tuple(r[k]) for r in rules["reactions"] for k in ("source", "target")})
