"""Build and analyze a reaction network from energy tables and rules."""
from dataclasses import asdict
from pathlib import Path
from .analysis import best_path, maximin_path, minimax_path, path_statistics, prune_network
from .data import load_energies
from .io import read_config, read_csv, resolve_path, sha256, write_csv, write_json
from .network import build_network
from .plotting import export_dot


def run_network(path):
    config = read_config(path)
    base = Path(path).resolve().parent
    root = resolve_path(base, config["output"])
    options = dict(config.get("network", {}))
    catalyst, facet = options.pop("catalyst", ""), options.pop("facet", "")
    if isinstance(config["energies"], dict):
        energies = dict(config["energies"])
    else:
        records, _ = load_energies(resolve_path(base, config["energies"]), catalyst, facet)
        energies = {key: row["free_energy"] for key, row in records.items()}
    if config.get("barriers"):
        for row in read_csv(resolve_path(base, config["barriers"]), ("catalyst", "source", "target", "barrier")):
            if row["catalyst"] != catalyst or row.get("facet", "") != facet:
                continue
            for key in ("potential", "ph", "temperature"):
                if row.get(key) and float(row[key]) != options.get(key, {"potential": 0, "ph": 0, "temperature": 298.15}[key]):
                    raise ValueError("Barrier condition mismatch: " + key)
            if row.get("potential_scale") and row["potential_scale"] != options.get("potential_scale", "RHE"):
                raise ValueError("Barrier potential scale mismatch")
            reaction = row["source"] + " -> " + row["target"]
            if reaction in energies:
                raise ValueError("Duplicate reaction entry: " + reaction)
            energies[reaction] = {"barrier": float(row["barrier"])}
    if "rules" in config:
        options["rules"] = config["rules"]
    elif config.get("rules_path"):
        options["rules"] = read_config(resolve_path(base, config["rules_path"]))
    graph = build_network(energies, metadata={"catalyst": catalyst, "facet": facet}, **options)
    root.mkdir(parents=True, exist_ok=True)
    graph.save(root / "network.json")
    export_dot(graph, root / "network.dot")
    write_csv(root / "states.csv", [dict(asdict(s), species="+".join(s.species)) for s in graph.states],
              ["id", "species", "free_energy", "product"])
    write_csv(root / "reactions.csv", [asdict(r) for r in graph.reactions],
              ["id", "source", "target", "kind", "delta_g", "barrier"])
    analysis = config.get("analysis", {})
    report = {"states": len(graph.states), "reactions": len(graph.reactions),
              "input_config_sha256": sha256(path), "metadata": graph.metadata, "products": {}}
    source = analysis.get("source")
    for i, target in enumerate(analysis.get("targets", [])):
        if analysis.get("weight"):
            objective = analysis.get("objective", "minimax")
            if objective not in ("minimax", "maximin"):
                raise ValueError("Raw weight search supports minimax or maximin")
            search = maximin_path if objective == "maximin" else minimax_path
            result = search(graph, source, target, weights=analysis["weight"])
        else:
            result = best_path(graph, source, target, mode=analysis.get("mode", "thermodynamic"),
                               objective=analysis.get("objective", "minimax"))
        statistics = path_statistics(graph, source, target, analysis.get("max_depth", 30),
                                     analysis.get("max_paths", 100000)) if analysis.get("path_counts", True) else None
        report["products"][target] = {"best_path": result, "path_statistics": statistics}
        if result["reachable"]:
            productive = prune_network(graph, source, [target], analysis.get("threshold"),
                                       mode=analysis.get("mode", "thermodynamic"))
            productive.save(root / ("product_{:03d}.json".format(i)))
            if config.get("plots", False):
                from .plotting import energy_diagram, network_plot
                nodes = {s.id: s for s in graph.states}
                if all(nodes[sid].free_energy is not None for sid in result["states"]):
                    energy_diagram(graph, result, root / ("energy_{:03d}.svg".format(i)))
                network_plot(productive, root / ("network_{:03d}.svg".format(i)),
                             (statistics or {}).get("node_counts"))
    if analysis.get("kinetics"):
        from .kinetics import absorption_probabilities, propagate, steady_state
        kinetic = analysis["kinetics"]
        kwargs = {k: kinetic[k] for k in ("temperature", "prefactor") if k in kinetic}
        report["kinetics"] = absorption_probabilities(graph, source, kinetic["products"], **kwargs)
        if kinetic.get("steady_state"):
            report["steady_state"] = steady_state(graph, **kwargs)
        if kinetic.get("times"):
            report["time_evolution"] = propagate(graph, {source: 1.0}, kinetic["times"],
                                                 products=kinetic["products"], **kwargs)
    write_json(root / "report.json", report)
    return report
