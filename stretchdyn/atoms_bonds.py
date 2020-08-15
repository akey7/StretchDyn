from dataclasses import dataclass, field
from typing import List, Dict
import numpy as np  # type: ignore
from numpy.linalg import norm  # type: ignore


class Atom:
    """
    This empty class is so I can get type checking to work down below.
    """
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
    
    @property
    def stretch_derivative(self):
        return self.k_ab * (2 * norm(self.r_ab) - 2 * self.r_e_ab) / 2

    @property
    def unit(self):
        return self.r_ab / norm(self.r_ab)

    @property
    def stretch_force(self):
        return -self.stretch_derivative * self.unit


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
