"""Generate synthetic energies covering all 629 legacy states; NOT research data."""
import csv
from pathlib import Path
from co2rr_network.chemistry import composition, electron_count
from co2rr_network.rules import enumerate_states, load_rules


def main():
    species = sorted({p for state in enumerate_states(load_rules("direct")) for p in state})
    root = Path(__file__).resolve().parent
    with (root / "synthetic_energies.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["catalyst", "facet", "species", "free_energy", "structure_id", "converged"])
        for i, name in enumerate(species):
            n = composition(name)
            # Deterministic artificial values used only to test graph coverage.
            energy = -.04 * electron_count(n) + .025 * n.get("O", 0) + (i % 7) * .01
            writer.writerow(["DEMO", "111", name, round(energy, 6), "synthetic_{}".format(i), "true"])
    print("Wrote synthetic energies for {} species".format(len(species)))


if __name__ == "__main__":
    main()
