"""Component recognition and dependency-free COO topology construction."""
from typing import List, Set, Tuple, Union
from .rules import load_rules
from .schema import state_id

Components = Union[List[str], Set[str]]


def _component_pool(values):
    if isinstance(values, str):
        raise ValueError("Components must be a list or set of single-component strings")
    pool = set()
    for value in values:
        if not isinstance(value, str):
            raise ValueError("Each component must be a string")
        component = value.strip()
        if not component or any(token in component for token in ("+", "->", "→")):
            raise ValueError("Expected a single component, got {!r}".format(value))
        pool.add(component)
    return pool


def recognize_reaction_pairs(comp_src: Components, comp_dst: Components, *,
                             profile="direct", rules=None) -> Set[str]:
    """Recognize allowed unit reactions between two component inventories.

    Each input is a pool of available single-component names, not the
    stoichiometric sides of one reaction. Thus {'CO'} can supply CO+CO.
    Repetition is retained in the recognized reaction, while duplicate input
    names do not change the pool. Ordered sites from the rules are preserved.
    Unknown components produce no matches; no reaction is inferred from formula
    balance alone. Custom rules may extend the built-in CO2RR mechanism.
    """
    source, target = _component_pool(comp_src), _component_pool(comp_dst)
    if not source or not target:
        return set()
    catalog = load_rules(profile, data=rules) if rules is not None else load_rules(profile)
    return {state_id(r["source"]) + " -> " + state_id(r["target"])
            for r in catalog["reactions"]
            if all(c in source for c in r["source"]) and all(c in target for c in r["target"])}


def _reaction_endpoints(pair):
    if not isinstance(pair, str):
        raise ValueError("Reaction pairs must be strings")
    pair = pair.replace("→", "->")
    if pair.count("->") != 1:
        raise ValueError("Expected 'source -> target', got {!r}".format(pair))
    sides = []
    for side in pair.split("->"):
        sides.append(state_id(tuple(part.strip() for part in side.split("+"))))
    return tuple(sides)


def make_reaction_graph(pairs: Set[str], *, return_labels=False):
    """Return (indices,) for a directed COO graph with logical shape (2, E).

    indices is [source_indices, target_indices], two Python integer lists.
    Nodes are whole reaction-side states (CO+CO is one node). Node labels and
    edges are sorted lexicographically so set iteration never affects indexing.
    With return_labels=True, return (indices, node_labels), where node_labels[i]
    decodes node i. Normalized duplicate edges collapse, as required by Set[str].
    No reverse edges or energy attributes are added. Empty input gives ([[], []],).
    """
    if isinstance(pairs, str):
        raise ValueError("pairs must be a collection of reaction strings")
    edges = sorted({_reaction_endpoints(pair) for pair in pairs})
    labels = tuple(sorted({node for edge in edges for node in edge}))
    node_index = {label: i for i, label in enumerate(labels)}
    indices = [[node_index[source] for source, _ in edges],
               [node_index[target] for _, target in edges]]
    return (indices, labels) if return_labels else (indices,)
