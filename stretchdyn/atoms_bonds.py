from dataclasses import dataclass, field
from typing import Any, Dict, List
import numpy as np  # type: ignore
from numpy.linalg import norm  # type: ignore


@dataclass
class Atom:
    """
    Represents a single atom in a molecule

    Instance attributes
    -------------------
    symbol: str
        The element symbol of the atom ("H", "Cl", etc)

    mass_amu: float
        The mass in atomic mass units.

    pos: np.array
        The position vector of the atom in angstroms.

    vel: np.array
        The velocity vector of the atom.

    bonds: Dict[str, Bond]
        The bonds in the into this atom. Note that the key is not the same
        as the element symbol. It needs to be unique. Note that below it
        is defined as type Dict[str, Any] because Bond is defined below.

        In fact there are a number of uses of bond that appear in this class.
        But, since Bond is defined below this class, they cannot be type hinted
        as Bond.
    """

    symbol: str
    mass_amu: float
    pos: np.array
    vel: np.array
    prev_accel: np.array
    bonds: Dict[str, Any] = field(default_factory=dict)
    pos_history: List[np.array] = field(default_factory=list)

    @property
    def net_stretch_force(self) -> np.array:
        """
        Unit returned: amu * Å / fs^2

        Returns
        -------
        np.array
            Net stretch force.
        """
        total_stretch_force = np.zeros(3)
        for bond in self.bonds.values():
            total_stretch_force -= bond.stretch_force
        return total_stretch_force

    @property
    def stretch_accel(self) -> np.array:
        """
        Units returned Å/fs^2

        Returns
        -------
        np.array
            The net acceleration due to stretch forces.
        """
        return self.net_stretch_force / self.mass_amu

    def update_pos_vel(self) -> None:
        """
        Using the acceleration, computes the new velocity and position.
        Uses velocity verlet integration.

        It doesn't return anything. Rather, it mutates the instnace
        variables.

        Returns
        -------
        None
            Simply mutates the instance variables.
        """
        dt_fs = 1
        next_accel = self.stretch_accel
        v_half_delta_t = self.vel + 0.5 * self.prev_accel * dt_fs
        self.pos = self.pos + v_half_delta_t * dt_fs
        self.vel = v_half_delta_t + 0.5 * next_accel * dt_fs
        self.prev_accel = next_accel
        self.pos_history.append(np.copy(self.pos))


@dataclass
class Bond:
    """
    This class represents a bond. It is comprised of instance attributes
    and properties that compute the properties of the bond needed for
    the molecular dynamics algorithm.

    Instance Attributes
    -------------------
    atom_a: Atom
        The first atom (atom A) in the bond.

    atom_b: Atom
        The second atom (atom B) in the bond

    r_e_ab: float
        The equilibrium distance between atoms A and B.

    k_ab: float
        The force constant for the stretch energy in the bond.
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
class Molecule:
    """
    This class just holds atoms, and has a convenience method
    update_all_atoms() which calculates the next positions
    and velocities for all the atoms.
    """

    atoms: Dict[str, Atom] = field(default_factory=dict)

    def update_all_atoms(self) -> None:
        """
        Returns
        -------
        None
            Returns nothing. It just mutates the atoms in place.
        """
        for atom in self.atoms.values():
            atom.update_pos_vel()
