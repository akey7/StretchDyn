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
    """
    This class represents a bond. It is comprised of instance attributes
    and properties that compute the properties of the bond needed for
    the molecular dynamics algorithm.
    """
    atom_a: Atom
    atom_b: Atom
    r_e_ab: float
    k_ab: float

    @property
    def r_ab(self) -> np.array:
        """
        The difference in positions between atoms B and A as a vector.

        Returns
        -------
        np.array
            Difference as a NumPy array
        """
        return self.atom_a.pos - self.atom_b.pos
    
    @property
    def stretch_derivative(self) -> float:
        """
        The derivative of the stretch energy over the distance between atoms
        A and B.

        Returns
        -------
        float
            The derivative of the stretch energy.
        """
        return self.k_ab * (2 * norm(self.r_ab) - 2 * self.r_e_ab) / 2

    @property
    def unit(self) -> np.array:
        """
        The unit vector in the direction from A to B.

        Returns
        -------
        np.array
            The unit vector in the direction pointing from A to B.
        """
        return self.r_ab / norm(self.r_ab)

    @property
    def stretch_force(self) -> np.array:
        """
        The force exerted on A by B.

        Returns
        -------
        np.array
            The force vector
        """
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
