"""Formula bookkeeping; names encode fragments, never molecular connectivity."""
import math
import re
from collections import Counter

KB_EV = 8.617333262145e-5
H_EV_S = 4.135667696e-15


def composition(name):
    """Count CHO including legacy O_CH4 / OH_CH3CH3 composite names."""
    name = name[:-2] if name.endswith("_d") else name
    result = Counter()
    for fragment in name.split("_"):
        tokens = re.findall(r"([A-Z][a-z]?)([0-9]*)", fragment)
        if not tokens or "".join(e + n for e, n in tokens) != fragment:
            raise ValueError("Invalid formula {!r}; use an explicit formula mapping".format(name))
        for element, number in tokens:
            if element not in ("C", "H", "O"):
                raise ValueError("CHO bookkeeping does not support {}".format(element))
            count = int(number or 1)
            if count <= 0:
                raise ValueError("Formula counts must be positive")
            result[element] += count
    return dict(result)


def state_composition(parts, formulas=None):
    result = Counter()
    formulas = formulas or {}
    for part in parts:
        result.update(composition(formulas.get(part, part)))
    return dict(result)


def electron_count(counts):
    return 4 * counts.get("C", 0) + counts.get("H", 0) - 2 * counts.get("O", 0)


def carbon_valence(name):
    counts = composition(name)
    if not counts.get("C"):
        raise ValueError("Carbon oxidation state needs at least one C atom")
    return (2 * counts.get("O", 0) - counts.get("H", 0)) / counts["C"]


def che_shift(counts, potential=0.0, ph=0.0, temperature=298.15, scale="RHE"):
    if temperature <= 0 or not all(math.isfinite(x) for x in (potential, ph, temperature)):
        raise ValueError("Finite conditions and temperature > 0 are required")
    if scale not in ("SHE", "RHE"):
        raise ValueError("potential_scale must be SHE or RHE")
    u = potential + (KB_EV * temperature * math.log(10) * ph if scale == "SHE" else 0.0)
    return electron_count(counts) * u


def balanced_step(source, target, formulas=None):
    left = state_composition(source, formulas)
    right = state_composition(target, formulas)
    if left.get("C", 0) != right.get("C", 0):
        raise ValueError("Carbon is not conserved: {} -> {}".format(source, target))
    water = left.get("O", 0) - right.get("O", 0)
    electrons = right.get("H", 0) - left.get("H", 0) + 2 * water
    return {"electrons_consumed": electrons, "protons_consumed": electrons,
            "water_produced": water}
