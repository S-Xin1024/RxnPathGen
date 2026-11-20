import math
import random
import unittest
from co2rr_network import Reaction, ReactionNetwork, State, maximin_path, minimax_path


def graph(edges):
    ids = sorted({node for source, target, _ in edges for node in (source, target)})
    network = ReactionNetwork([State(s, (s,)) for s in ids],
        [Reaction(str(i), a, b, "edge") for i, (a, b, _) in enumerate(edges)])
    return network, {str(i): w for i, (_, _, w) in enumerate(edges)}


class BottleneckTests(unittest.TestCase):
    def test_minimax_is_not_shortest_sum(self):
        g, w = graph([("C", "CO", 6), ("CO", "CO2", 6), ("C", "CO2", 10)])
        result = minimax_path(g, "C", "CO2", w)
        self.assertEqual(result["cost"], 6)
        self.assertEqual(result["states"], ["C", "CO", "CO2"])

    def test_minimax_preserves_negative_weights(self):
        g, w = graph([("C", "CO", -5), ("CO", "CO2", -3), ("C", "CO2", -1)])
        self.assertEqual(minimax_path(g, "C", "CO2", w)["cost"], -3)

    def test_maximin_and_parallel_channels(self):
        g, w = graph([("C", "CO", 2), ("C", "CO", 8), ("CO", "CO2", 7), ("C", "CO2", 6)])
        self.assertEqual(maximin_path(g, "C", "CO2", w)["cost"], 7)
        self.assertEqual(maximin_path(g, "C", "CO2", w)["reactions"], ["1", "2"])

    def test_empty_and_unreachable(self):
        g, w = graph([("C", "CO", -2)])
        self.assertIsNone(minimax_path(g, "C", "C", w)["cost"])
        self.assertTrue(minimax_path(g, "C", "C", w)["empty_path"])
        self.assertFalse(maximin_path(g, "CO", "C", w)["reachable"])

    def test_callable_and_invalid_weights(self):
        g, w = graph([("C", "CO", 3)])
        self.assertEqual(minimax_path(g, "C", "CO", lambda r: w[r.id])["cost"], 3)
        for values in ({}, {"0": float("nan")}, {"0": float("inf")}):
            with self.assertRaises(ValueError):
                minimax_path(g, "C", "CO", values)

    def test_matches_exhaustive_simple_paths_on_small_cyclic_graphs(self):
        rng = random.Random(42)
        for _ in range(40):
            nodes = ["C", "CO", "CO2", "CH", "CH2"]
            edges = [(a, b, rng.randint(-8, 8)) for a in nodes for b in nodes
                     if a != b and rng.random() < .35]
            edges += [("C", "CO2", rng.randint(-8, 8))]
            g, w = graph(edges)
            out = {s.id: [] for s in g.states}
            for r in g.reactions:
                out[r.source].append(r)
            costs = []
            stack = [("C", {"C"}, [])]
            while stack:
                node, seen, weights = stack.pop()
                if node == "CO2":
                    costs.append(weights)
                    continue
                for r in out[node]:
                    if r.target not in seen:
                        stack.append((r.target, seen | {r.target}, weights + [w[r.id]]))
            self.assertEqual(minimax_path(g, "C", "CO2", w)["cost"], min(max(v) for v in costs))
            self.assertEqual(maximin_path(g, "C", "CO2", w)["cost"], max(min(v) for v in costs))
