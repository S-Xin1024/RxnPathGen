"""Commands for graph construction, analysis and path search."""
import argparse
import json
import sys
from .io import write_json


def parser():
    p = argparse.ArgumentParser(prog="co2rr", description="Reaction network construction and analysis")
    p.add_argument("--version", action="version", version="co2rr-network 0.2.0")
    sub = p.add_subparsers(dest="command")
    network = sub.add_parser("network", help="Build and analyze a network from JSON configuration")
    network.add_argument("config")
    rules = sub.add_parser("rules", help="Export or summarize the reaction rules")
    rules.add_argument("--profile", choices=("direct", "out"), default="direct")
    rules.add_argument("--output")
    legacy = sub.add_parser("import-legacy", help="Convert catalyst_species,value energy CSV")
    legacy.add_argument("input")
    legacy.add_argument("output")
    legacy.add_argument("--facet", default="")
    path = sub.add_parser("path", help="Search an existing network")
    path.add_argument("network")
    path.add_argument("source")
    path.add_argument("target")
    path.add_argument("--mode", choices=("thermodynamic", "kinetic"), default="thermodynamic")
    path.add_argument("--objective", choices=("minimax", "sum", "maximin"), default="minimax")
    path.add_argument("--weight", choices=("delta_g", "barrier"),
                      help="Use raw edge weights; minimax and maximin accept negative values")
    stats = sub.add_parser("paths", help="Count paths and node occurrences")
    stats.add_argument("network")
    stats.add_argument("source")
    stats.add_argument("target")
    stats.add_argument("--max-depth", type=int, default=30)
    stats.add_argument("--max-paths", type=int, default=100000)
    prune = sub.add_parser("prune", help="Export the source-to-target productive subnetwork")
    prune.add_argument("network")
    prune.add_argument("source")
    prune.add_argument("targets", nargs="+")
    prune.add_argument("--output", required=True)
    prune.add_argument("--threshold", type=float)
    prune.add_argument("--mode", choices=("thermodynamic", "kinetic"), default="thermodynamic")
    plot = sub.add_parser("plot", help="Plot a network or energy table")
    plot.add_argument("input")
    plot.add_argument("output")
    plot.add_argument("--kind", choices=("network", "heatmap"), default="network")
    plot.add_argument("--facet", default="")
    return p


def execute(args):
    from .schema import ReactionNetwork
    cmd = args.command
    if cmd == "network":
        from .workflow import run_network
        result = run_network(args.config)
        return {"states": result["states"], "reactions": result["reactions"],
                "products": {k: v["best_path"] for k, v in result["products"].items()}}
    if cmd == "rules":
        from .rules import enumerate_states, load_rules
        rules = load_rules(args.profile)
        if args.output:
            write_json(args.output, rules)
        return {"profile": args.profile, "states": len(enumerate_states(rules)),
                "rules": len(rules["reactions"])}
    if cmd == "import-legacy":
        from .data import import_legacy_csv
        return {"rows": import_legacy_csv(args.input, args.output, facet=args.facet)}
    if cmd == "path":
        from .analysis import best_path, maximin_path, minimax_path
        network = ReactionNetwork.load(args.network)
        if args.weight:
            if args.objective == "sum":
                raise ValueError("Raw signed weights are supported for minimax/maximin, not sum")
            search = maximin_path if args.objective == "maximin" else minimax_path
            return search(network, args.source, args.target, weights=args.weight)
        if args.objective == "maximin":
            raise ValueError("maximin requires an explicit --weight")
        return best_path(network, args.source, args.target, args.mode, args.objective)
    if cmd == "paths":
        from .analysis import path_statistics
        return path_statistics(ReactionNetwork.load(args.network), args.source, args.target,
                               args.max_depth, args.max_paths)
    if cmd == "prune":
        from .analysis import prune_network
        graph = prune_network(ReactionNetwork.load(args.network), args.source, args.targets,
                              args.threshold, args.mode)
        graph.save(args.output)
        return {"output": args.output, "states": len(graph.states), "reactions": len(graph.reactions)}
    if cmd == "plot":
        from .plotting import heatmap, network_plot
        if args.kind == "heatmap":
            heatmap(args.input, args.output, args.facet)
        else:
            network_plot(ReactionNetwork.load(args.input), args.output)
        return {"output": args.output}
    raise ValueError("Unknown command")


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    if args.command is None:
        p.print_help()
        return 0
    try:
        result = execute(args)
    except (ValueError, OSError, RuntimeError, KeyError, ImportError) as exc:
        print("co2rr: {}".format(exc), file=sys.stderr)
        return 2
    print(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False))
    return 0
