import importlib.util
import math
import unittest
from co2rr_network.schema import State, Reaction, ReactionNetwork

NUMPY = importlib.util.find_spec("numpy") is not None
SCIPY = importlib.util.find_spec("scipy") is not None


def graph():
    return ReactionNetwork([State("CO+CO", ("CO", "CO"), 0), State("COCO", ("COCO",), -.1),
                            State("CHOCHO_d", ("CHOCHO_d",), -.2, True)],
                           [Reaction("a", "CO+CO", "COCO", "cc", -.1, .3),
                            Reaction("b", "COCO", "CO+CO", "cc", .1, .4),
                            Reaction("c", "COCO", "CHOCHO_d", "h", -.1, .3)])


@unittest.skipUnless(NUMPY, "optional numpy unavailable")
class KineticsTests(unittest.TestCase):
    def test_probability_and_rate_conservation(self):
        import numpy as np
        from co2rr_network.kinetics import generator, jump_matrix
        _, p = jump_matrix(graph(), ["CHOCHO_d"])
        _, q = generator(graph(), ["CHOCHO_d"])
        np.testing.assert_allclose(p.sum(axis=1), 1)
        np.testing.assert_allclose(q.sum(axis=1), 0, atol=1e-8)

    def test_nonproductive_trap(self):
        from co2rr_network.kinetics import absorption_probabilities
        g = graph()
        g.states.append(State("CHO", ("CHO",), -.3))
        g.reactions = [Reaction("prod", "CO+CO", "CHOCHO_d", "h", -.2, .3),
                       Reaction("trap", "CO+CO", "CHO", "h", -.3, .3)]
        result = absorption_probabilities(g, "CO+CO", ["CHOCHO_d"])
        self.assertAlmostEqual(result["trapped_probability"], .5)

    def test_detailed_balance(self):
        from co2rr_network.kinetics import steady_state
        from co2rr_network.chemistry import KB_EV
        g = graph()
        g.states, g.reactions = g.states[:2], g.reactions[:2]
        p = steady_state(g)
        self.assertAlmostEqual(p["COCO"] / p["CO+CO"], math.exp(.1 / (KB_EV * 298.15)), places=7)

    def test_nonunique_steady_state(self):
        from co2rr_network.kinetics import steady_state
        g = graph()
        g.reactions = []
        with self.assertRaises(ValueError):
            steady_state(g)

    def test_high_barriers_and_unknown_reaction_energy(self):
        from co2rr_network import build_network
        from co2rr_network.kinetics import jump_matrix
        g = build_network({"CO+CO -> COCO": {"barrier": 100.5}})
        _, p = jump_matrix(g)
        self.assertAlmostEqual(p[0].sum(), 1)

    @unittest.skipUnless(SCIPY, "optional scipy unavailable")
    def test_time_evolution_conserves_probability(self):
        from co2rr_network.kinetics import propagate
        result = propagate(graph(), {"CO+CO": 1}, [0, 1e-7, 1e-5], products=["CHOCHO_d"])
        for row in result["probabilities"]:
            self.assertAlmostEqual(sum(row), 1)
