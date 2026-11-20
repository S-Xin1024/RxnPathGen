"""Dependency-free path search, pruning and bounded path statistics."""
import heapq
import math
from collections import Counter, defaultdict, deque
from .schema import ReactionNetwork
from .io import finite


def adjacency(network):
    result = {s.id: [] for s in network.states}
    for reaction in network.reactions:
        result[reaction.source].append(reaction)
    return result


def edge_cost(reaction, mode):
    if mode == "thermodynamic":
        if reaction.delta_g is None:
            raise ValueError("Thermodynamic analysis requires delta_g for {}".format(reaction.id))
        return max(0.0, reaction.delta_g)
    if mode == "kinetic":
        if reaction.barrier is None:
            raise ValueError("Kinetic analysis requires a barrier for {} ({} -> {})".format(
                reaction.id, reaction.source, reaction.target))
        return reaction.barrier
    raise ValueError("mode must be thermodynamic or kinetic")


def _bottleneck_path(network, source, target, weights, maximize):
    adj = adjacency(network)
    if source not in adj or target not in adj:
        raise ValueError("Unknown source or target")
    values = {}
    for r in network.reactions:
        if isinstance(weights, str):
            if weights not in ("delta_g", "barrier"):
                raise ValueError("Named weights must be delta_g or barrier")
            value = getattr(r, weights)
        elif callable(weights):
            value = weights(r)
        else:
            if r.id not in weights:
                raise ValueError("Missing weight for {}".format(r.id))
            value = weights[r.id]
        values[r.id] = finite(value, "weight of " + r.id)
    objective = "maximin" if maximize else "minimax"
    result = {"source": source, "target": target, "objective": objective,
              "weight": weights if isinstance(weights, str) else "custom"}
    if source == target:
        return dict(result, reachable=True, empty_path=True, cost=None,
                    states=[source], reactions=[], step_costs=[], bottleneck_reactions=[])
    identity = math.inf if maximize else -math.inf
    scores, previous = {source: identity}, {}
    queue = [(-identity if maximize else identity, source)]
    while queue:
        priority, node = heapq.heappop(queue)
        score = -priority if maximize else priority
        if score != scores[node]:
            continue
        if node == target:
            break
        for r in adj[node]:
            candidate = min(score, values[r.id]) if maximize else max(score, values[r.id])
            default = -math.inf if maximize else math.inf
            better = candidate > scores.get(r.target, default) if maximize else candidate < scores.get(r.target, default)
            if better:
                scores[r.target] = candidate
                previous[r.target] = r
                heapq.heappush(queue, (-candidate if maximize else candidate, r.target))
    if target not in scores:
        return dict(result, reachable=False, empty_path=False, cost=None,
                    states=[], reactions=[], step_costs=[], bottleneck_reactions=[])
    nodes, steps, node = [target], [], target
    while node != source:
        r = previous[node]
        steps.append(r)
        node = r.source
        nodes.append(node)
    nodes.reverse()
    steps.reverse()
    costs = [values[r.id] for r in steps]
    return dict(result, reachable=True, empty_path=False, cost=scores[target], states=nodes,
                reactions=[r.id for r in steps], step_costs=costs,
                bottleneck_reactions=[r.id for r in steps if abs(values[r.id] - scores[target]) < 1e-9])


def minimax_path(network, source, target, weights="delta_g"):
    """Minimize the largest raw edge weight, including finite negative weights.

    weights is an edge-attribute name, a mapping from reaction ID to weight,
    or a callable receiving a Reaction. No clipping or energy shift is applied.
    An empty path has cost=None because its maximum edge weight is undefined.
    """
    return _bottleneck_path(network, source, target, weights, maximize=False)


def maximin_path(network, source, target, weights="delta_g"):
    """Maximize the smallest raw edge weight (widest/bottleneck-capacity path)."""
    return _bottleneck_path(network, source, target, weights, maximize=True)


def best_path(network, source, target, mode="thermodynamic", objective="minimax"):
    """Minimize maximum step cost or sum of nonnegative step costs.

    Thermodynamic costs use max(0,dG), never an arbitrary constant shift.
    This is a path-screening descriptor, not a kinetic rate calculation.
    """
    if objective not in ("minimax", "sum"):
        raise ValueError("objective must be minimax or sum")
    adj = adjacency(network)
    if source not in adj or target not in adj:
        raise ValueError("Unknown source or target")
    distances = {source: 0.0}
    previous, queue = {}, [(0.0, source)]
    while queue:
        distance, node = heapq.heappop(queue)
        if distance != distances[node]:
            continue
        if node == target:
            break
        for reaction in adj[node]:
            weight = edge_cost(reaction, mode)
            candidate = max(distance, weight) if objective == "minimax" else distance + weight
            if candidate < distances.get(reaction.target, math.inf):
                distances[reaction.target] = candidate
                previous[reaction.target] = reaction
                heapq.heappush(queue, (candidate, reaction.target))
    if target not in distances:
        return {"reachable": False, "source": source, "target": target, "cost": None,
                "states": [], "reactions": [], "mode": mode, "objective": objective}
    nodes, steps, current = [target], [], target
    while current != source:
        reaction = previous[current]
        steps.append(reaction)
        current = reaction.source
        nodes.append(current)
    nodes.reverse()
    steps.reverse()
    costs = [edge_cost(r, mode) for r in steps]
    bottleneck = max(costs, default=0.0)
    return {"reachable": True, "source": source, "target": target, "cost": distances[target],
            "states": nodes, "reactions": [r.id for r in steps], "mode": mode,
            "objective": objective, "step_costs": costs,
            "bottleneck_reactions": [r.id for r, c in zip(steps, costs) if abs(c - bottleneck) < 1e-9],
            "interpretation": "path descriptor; does not identify a kinetic rate-determining step"}


def reachable(adj, starts):
    seen, stack = set(starts), list(starts)
    while stack:
        for target in adj.get(stack.pop(), []):
            if target not in seen:
                seen.add(target)
                stack.append(target)
    return seen


def prune_network(network, source, targets, threshold=None, mode="thermodynamic"):
    nodes = {s.id for s in network.states}
    if source not in nodes or not set(targets).issubset(nodes):
        raise ValueError("Unknown source or targets")
    edges = [r for r in network.reactions if threshold is None or edge_cost(r, mode) <= threshold]
    out, incoming = defaultdict(list), defaultdict(list)
    for r in edges:
        out[r.source].append(r.target)
        incoming[r.target].append(r.source)
    keep = reachable(out, [source]) & reachable(incoming, targets)
    return ReactionNetwork([s for s in network.states if s.id in keep],
                           [r for r in edges if r.source in keep and r.target in keep],
                           dict(network.metadata, pruning={"source": source, "targets": list(targets),
                                                         "threshold": threshold, "mode": mode}))


def path_statistics(network, source, target, max_depth=30, max_paths=100000):
    """Edge-distinct paths and node occurrence. DAGs use exact integer DP.

    On cyclic graphs only simple paths up to max_depth edges are considered,
    with max_paths as a work bound. Truncation is reported, never hidden.
    """
    if max_depth < 0 or max_paths < 1:
        raise ValueError("max_depth >= 0 and max_paths >= 1 required")
    graph = prune_network(network, source, [target])
    if not graph.states:
        return {"total_paths": 0, "node_counts": {}, "exact": True, "method": "unreachable"}
    adj = adjacency(graph)
    # A path terminates on reaching target, even if target has outgoing edges.
    adj[target] = []
    indegree = {s.id: 0 for s in graph.states}
    for values in adj.values():
        for r in values:
            indegree[r.target] += 1
    queue = deque(sorted(n for n, degree in indegree.items() if not degree))
    order = []
    while queue:
        node = queue.popleft()
        order.append(node)
        for r in adj[node]:
            indegree[r.target] -= 1
            if not indegree[r.target]:
                queue.append(r.target)
    if len(order) == len(adj):
        forward, backward = defaultdict(int), defaultdict(int)
        forward[source], backward[target] = 1, 1
        for node in order:
            for r in adj[node]:
                forward[r.target] += forward[node]
        for node in reversed(order):
            for r in adj[node]:
                backward[node] += backward[r.target]
        return {"total_paths": forward[target], "node_counts": {
            n: forward[n] * backward[n] for n in order}, "exact": True, "method": "dag_dynamic_programming"}
    counts, total, depth_limited, count_limited = Counter(), 0, False, False
    expansions, expansion_limited = 0, False
    max_expansions = max_paths * max(1, max_depth)
    stack = [(source, (source,))]
    while stack:
        if expansions >= max_expansions:
            expansion_limited = True
            break
        expansions += 1
        node, path = stack.pop()
        if node == target:
            if total == max_paths:
                count_limited = True
                break
            total += 1
            counts.update(path)
            continue
        choices = [r for r in adj[node] if r.target not in path]
        if len(path) - 1 >= max_depth:
            depth_limited |= bool(choices)
            continue
        for reaction in reversed(choices):
            stack.append((reaction.target, path + (reaction.target,)))
    return {"total_paths": total, "node_counts": dict(counts),
            "exact": not (depth_limited or count_limited or expansion_limited), "depth_limited": depth_limited,
            "count_limited": count_limited, "method": "bounded_simple_paths",
            "expansion_limited": expansion_limited, "expansions": expansions,
            "max_depth": max_depth, "max_paths": max_paths}
