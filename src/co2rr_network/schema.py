"""Serializable state graph. A node represents a whole surface configuration."""
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Tuple
from .io import finite, read_config, write_json


def state_id(parts):
    if isinstance(parts, str):
        raise ValueError("species must be a sequence of names, not one string")
    parts = tuple(parts)
    if not parts or any(not isinstance(p, str) or not p or "+" in p for p in parts):
        raise ValueError("A state needs nonempty species names without '+'")
    return "+".join(parts)


@dataclass
class State:
    id: str
    species: Tuple[str, ...]
    free_energy: Optional[float] = None
    product: bool = False
    metadata: Dict = field(default_factory=dict)

    def __post_init__(self):
        state_id(self.species)
        self.species = tuple(self.species)
        if self.id != state_id(self.species):
            raise ValueError("State id must equal '+'.join(species)")
        if self.free_energy is not None:
            self.free_energy = finite(self.free_energy, "free_energy")


@dataclass
class Reaction:
    id: str
    source: str
    target: str
    kind: str
    delta_g: Optional[float] = None
    barrier: Optional[float] = None
    metadata: Dict = field(default_factory=dict)

    def __post_init__(self):
        if self.delta_g is not None:
            self.delta_g = finite(self.delta_g, "delta_g")
        if self.barrier is not None:
            self.barrier = finite(self.barrier, "barrier")
            if self.barrier < 0:
                raise ValueError("Activation barriers must be nonnegative")


@dataclass
class ReactionNetwork:
    states: List[State]
    reactions: List[Reaction]
    metadata: Dict = field(default_factory=dict)

    def __post_init__(self):
        nodes = {s.id for s in self.states}
        if len(nodes) != len(self.states):
            raise ValueError("Duplicate state IDs")
        if len({r.id for r in self.reactions}) != len(self.reactions):
            raise ValueError("Duplicate reaction IDs")
        energies = {s.id: s.free_energy for s in self.states}
        for reaction in self.reactions:
            if reaction.source not in nodes or reaction.target not in nodes:
                raise ValueError("Reaction references a missing state: {}".format(reaction.id))
            known = (energies[reaction.target], energies[reaction.source], reaction.delta_g)
            if all(value is not None for value in known) and abs(known[0] - known[1] - known[2]) > 1e-7:
                raise ValueError("delta_g inconsistent with node energies: {}".format(reaction.id))

    def save(self, path):
        write_json(path, {"schema_version": 1, "states": [asdict(s) for s in self.states],
                          "reactions": [asdict(r) for r in self.reactions], "metadata": self.metadata})

    @classmethod
    def load(cls, path):
        data = read_config(path)
        if data.get("schema_version") != 1:
            raise ValueError("Unsupported network schema")
        return cls([State(**s) for s in data["states"]],
                   [Reaction(**r) for r in data["reactions"]], data.get("metadata", {}))
