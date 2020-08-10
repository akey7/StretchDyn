from dataclasses import dataclass, field
from typing import List, Dict
import numpy as np  # type: ignore
from numpy.linalg import norm  # type: ignore


class Atom:
    pass


@dataclass
class Bond:
    atom_a: Atom
    atom_b: Atom
    r_e_ab: float
    k_ab: float

    @property
    def r_ab(self):
        return self.atom_a.pos - self.atom_b.pos


@dataclass
class Atom:
    symbol: str
    mass_amu: float
    pos: np.array
    vel: np.array
    bonds: Dict[str, Bond] = field(default_factory=dict)


@dataclass
class Molecule:
    atoms: Dict[str, Atom] = field(default_factory=dict)
