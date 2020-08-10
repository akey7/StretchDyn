from dataclasses import dataclass
import numpy as np


@dataclass
class Atom:
    mass_amu: float
    position: np.array


@dataclass
class Bond:
    atom_a: Atom
    atom_b: Atom
    r_e_ab: float
    k_ab: float

    @property
    def r_ab(self):
        return self.atom_a.position - self.atom_b.position
