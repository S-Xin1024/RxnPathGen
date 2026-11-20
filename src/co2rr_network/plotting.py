"""Exportable plots with lazy matplotlib import; no GPU/model dependencies."""
import json
from pathlib import Path
from .chemistry import electron_count


def _pyplot():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("Plotting requires: pip install 'co2rr-network[plots]'") from exc
    return plt


def energy_diagram(network, result, output):
    if not result.get("reachable"):
        raise ValueError("Cannot plot an unreachable path")
    plt = _pyplot()
    states = {s.id: s for s in network.states}
    reactions = {r.id: r for r in network.reactions}
    nodes = result["states"]
    if any(states[sid].free_energy is None for sid in nodes):
        raise ValueError("An absolute energy diagram needs state free_energy values")
    reference = states[nodes[0]].free_energy
    fig, ax = plt.subplots(figsize=(max(6, len(nodes) * 1.0), 4.5))
    for i, sid in enumerate(nodes):
        energy = states[sid].free_energy - reference
        ax.plot([i - .24, i + .24], [energy, energy], color="#2563a6", lw=2.5)
        if i:
            previous = states[nodes[i - 1]].free_energy - reference
            r = reactions[result["reactions"][i - 1]]
            if r.barrier is None:
                ax.plot([i - .76, i - .24], [previous, energy], "--", color="0.5", lw=1)
            else:
                ax.plot([i - .76, i - .5, i - .24], [previous, previous + r.barrier, energy],
                        color="#b24b32", lw=1.3)
    ax.set_xticks(range(len(nodes)))
    ax.set_xticklabels(nodes, rotation=40, ha="right")
    ax.set_ylabel("Relative formation free energy (eV)")
    ax.set_xlabel("Reaction coordinate")
    ax.axhline(0, color="0.8", lw=.6)
    fig.tight_layout()
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def network_plot(network, output, counts=None):
    plt = _pyplot()
    groups, positions = {}, {}
    for s in network.states:
        x = electron_count(s.metadata.get("composition", {}))
        groups.setdefault(x, []).append(s)
    for x, values in groups.items():
        for i, state in enumerate(sorted(values, key=lambda s: s.id)):
            positions[state.id] = (x, i - (len(values) - 1) / 2)
    fig, ax = plt.subplots(figsize=(14, max(5, min(18, len(network.states) / 20))))
    for r in network.reactions:
        x0, y0 = positions[r.source]
        x1, y1 = positions[r.target]
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops={"arrowstyle": "->", "color": "0.7", "lw": .4, "alpha": .6})
    counts = counts or {}
    import math
    has_energies = all(s.free_energy is not None for s in network.states)
    if counts:
        values = [math.log10(1 + counts.get(s.id, 0)) for s in network.states]
        label = "log10(1 + path count)"
    elif has_energies:
        values, label = [s.free_energy for s in network.states], "Formation G (eV)"
    else:
        values = [sum(r.source == s.id or r.target == s.id for r in network.reactions) for s in network.states]
        label = "Node degree"
    dots = ax.scatter([positions[s.id][0] for s in network.states], [positions[s.id][1] for s in network.states],
                      c=values, cmap="viridis", s=20, zorder=3)
    if len(network.states) <= 80:
        for s in network.states:
            ax.annotate(s.id, positions[s.id], xytext=(3, 3), textcoords="offset points", fontsize=6)
    fig.colorbar(dots, ax=ax, label=label)
    ax.set_xlabel("Electrons consumed relative to CO2/H2O")
    ax.set_ylabel("State index within electron count")
    fig.tight_layout()
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def heatmap(csv_path, output, facet=""):
    from .io import read_csv
    plt = _pyplot()
    import numpy as np
    rows = [r for r in read_csv(csv_path, ("catalyst", "species", "free_energy"))
            if r.get("facet", "") == facet and r.get("converged", "true").lower() not in ("false", "0", "no")]
    catalysts, species = sorted({r["catalyst"] for r in rows}), sorted({r["species"] for r in rows})
    if not rows:
        raise ValueError("No heatmap energies for selected facet")
    data = np.full((len(catalysts), len(species)), np.nan)
    for r in rows:
        i, j = catalysts.index(r["catalyst"]), species.index(r["species"])
        value = float(r["free_energy"])
        data[i, j] = min(data[i, j], value) if np.isfinite(data[i, j]) else value
    fig, ax = plt.subplots(figsize=(max(6, len(species) * .25), max(3, len(catalysts) * .35)))
    im = ax.imshow(np.ma.masked_invalid(data), aspect="auto", cmap="coolwarm")
    ax.set_xticks(range(len(species)))
    ax.set_xticklabels(species, rotation=90, fontsize=7)
    ax.set_yticks(range(len(catalysts)))
    ax.set_yticklabels(catalysts)
    fig.colorbar(im, ax=ax, label="Formation G at U=0, pH=0 (eV)")
    fig.tight_layout()
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def export_dot(network, output):
    lines = ["digraph CO2RR {", "rankdir=LR;"]
    for s in network.states:
        energy = "" if s.free_energy is None else "\\nG={:.4f} eV".format(s.free_energy)
        lines.append("{} [label={}, shape={}];".format(
            json.dumps(s.id), json.dumps(s.id + energy),
            "doublecircle" if s.product else "ellipse"))
    for r in network.reactions:
        label = r.kind
        if r.delta_g is not None:
            label += " dG={:.3f}".format(r.delta_g)
        if r.barrier is not None:
            label += " Ea={:.3f}".format(r.barrier)
        lines.append("{} -> {} [label={}];".format(json.dumps(r.source), json.dumps(r.target), json.dumps(label)))
    lines.append("}")
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text("\n".join(lines) + "\n", encoding="utf-8")
