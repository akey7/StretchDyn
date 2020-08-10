from dataclasses import dataclass  # type: ignore
import numpy as np  # type: ignore


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
    def r_ab(self) -> np.array:
        return self.atom_a.position - self.atom_b.position
