"""First-order state kinetics from transition-state-theory barriers.

This is a single-configuration CTMC, not a coverage-dependent microkinetic model.
Barrier and G differences must share conditions, reference and corrections.
"""
import math
from .chemistry import H_EV_S, KB_EV


def _numpy():
    try:
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("Kinetics requires: pip install 'co2rr-network[kinetics]'") from exc
    return np


def log_rates(network, temperature=None, prefactor=None):
    temperature = (network.metadata.get("conditions", {}).get("temperature", 298.15)
                   if temperature is None else temperature)
    if not math.isfinite(temperature) or temperature <= 0:
        raise ValueError("temperature must be finite and positive")
    prefactor = KB_EV * temperature / H_EV_S if prefactor is None else prefactor
    if not math.isfinite(prefactor) or prefactor <= 0:
        raise ValueError("prefactor must be finite and positive")
    logs = []
    for r in network.reactions:
        if r.barrier is None:
            raise ValueError("Missing barrier for {}; cannot infer kinetic rates from dG".format(r.id))
        if r.delta_g is not None and r.barrier + 1e-7 < max(0, r.delta_g):
            raise ValueError("Barrier below endpoint: {}".format(r.id))
        logs.append((r, math.log(prefactor) - r.barrier / (KB_EV * temperature)))
    return logs


def jump_matrix(network, products=(), **kwargs):
    """Embedded jump probabilities, calculated with row-wise log-sum-exp."""
    np = _numpy()
    ids = [s.id for s in network.states]
    index = {sid: i for i, sid in enumerate(ids)}
    if not set(products).issubset(index):
        raise ValueError("Unknown absorbing products")
    rows = {sid: [] for sid in ids}
    for r, logk in log_rates(network, **kwargs):
        if r.source != r.target:
            rows[r.source].append((r.target, logk))
    matrix = np.zeros((len(ids), len(ids)))
    for sid, values in rows.items():
        i = index[sid]
        if sid in products or not values:
            matrix[i, i] = 1.0
            continue
        shift = max(v for _, v in values)
        weights = [(target, math.exp(logk - shift)) for target, logk in values]
        norm = sum(v for _, v in weights)
        for target, value in weights:
            matrix[i, index[target]] += value / norm
    return ids, matrix


def generator(network, products=(), **kwargs):
    np = _numpy()
    ids = [s.id for s in network.states]
    index = {sid: i for i, sid in enumerate(ids)}
    if not set(products).issubset(index):
        raise ValueError("Unknown products")
    matrix = np.zeros((len(ids), len(ids)))
    for r, logk in log_rates(network, **kwargs):
        if r.source not in products and r.source != r.target:
            k = math.exp(logk)
            if k == 0:
                raise ValueError("Rate underflow for {}; use jump_matrix for absorption".format(r.id))
            matrix[index[r.source], index[r.target]] += k
    for i in range(len(ids)):
        matrix[i, i] = -matrix[i].sum()
    return ids, matrix


def absorption_probabilities(network, source, products, **kwargs):
    from .analysis import reachable
    np = _numpy()
    products = list(dict.fromkeys(products))
    if not products:
        raise ValueError("Supply absorbing products")
    ids, matrix = jump_matrix(network, products, **kwargs)
    if source not in ids:
        raise ValueError("Unknown source")
    index = {sid: i for i, sid in enumerate(ids)}
    reverse = {sid: [] for sid in ids}
    for i, sid in enumerate(ids):
        for j in np.flatnonzero(matrix[i]):
            if i != j:
                reverse[ids[j]].append(sid)
    productive = reachable(reverse, products)
    if source in products:
        probabilities = {p: float(p == source) for p in products}
    elif source not in productive:
        probabilities = dict.fromkeys(products, 0.0)
    else:
        transient = [index[sid] for sid in ids if sid in productive and sid not in products]
        absorbing = [index[sid] for sid in products]
        q = matrix[np.ix_(transient, transient)]
        r = matrix[np.ix_(transient, absorbing)]
        system = np.eye(len(q)) - q
        if np.linalg.cond(system) > 1e13:
            raise ValueError("Absorption system is ill-conditioned; rates span unresolved timescales")
        b = np.linalg.solve(system, r)
        values = b[transient.index(index[source])]
        if np.any(values < -1e-8) or values.sum() > 1 + 1e-8:
            raise ValueError("Numerically invalid absorption probabilities")
        probabilities = {p: float(max(0, min(1, v))) for p, v in zip(products, values)}
    return {"product_probabilities": probabilities,
            "trapped_probability": max(0.0, 1 - sum(probabilities.values())),
            "model": "first-order configuration CTMC; not Faradaic efficiency"}


def steady_state(network, **kwargs):
    """Unique stationary distribution of the continuous-time generator."""
    np = _numpy()
    ids, q = generator(network, **kwargs)
    scale = abs(q).max()
    if scale == 0 or np.linalg.matrix_rank(q / scale) != len(q) - 1:
        raise ValueError("Steady state is not uniquely numerically resolved (multiple closed classes, "
                         "disconnected sinks, or unresolved rate scales)")
    a = q.T / scale
    a[-1] = 1.0
    b = np.zeros(len(q))
    b[-1] = 1.0
    pi = np.linalg.solve(a, b)
    if pi.min() < -1e-8 or abs(pi.sum() - 1) > 1e-8:
        raise ValueError("Invalid stationary distribution")
    return {sid: float(max(0.0, v)) for sid, v in zip(ids, pi)}


def propagate(network, initial, times, products=(), **kwargs):
    np = _numpy()
    try:
        from scipy.linalg import expm
    except ImportError as exc:
        raise RuntimeError("Time propagation requires scipy") from exc
    ids, q = generator(network, products, **kwargs)
    if not set(initial).issubset(ids):
        raise ValueError("Unknown initial states")
    p0 = np.array([initial.get(sid, 0.0) for sid in ids], dtype=float)
    if not np.all(np.isfinite(p0)) or p0.min() < 0 or abs(p0.sum() - 1) > 1e-9:
        raise ValueError("Initial distribution must be finite, nonnegative and normalized")
    if any(not math.isfinite(t) or t < 0 for t in times):
        raise ValueError("Times must be finite and nonnegative (s)")
    values = [p0 @ expm(q * t) for t in times]
    if any(v.min() < -1e-8 or abs(v.sum() - 1) > 1e-7 for v in values):
        raise ValueError("Propagation lost probability conservation")
    return {"states": ids, "time_s": list(times), "probabilities": [v.tolist() for v in values]}
