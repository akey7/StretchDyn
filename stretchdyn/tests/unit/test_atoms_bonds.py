import pytest
import numpy as np
from stretchdyn.atoms_bonds import Atom, Bond


@pytest.fixture
def atom_cl():
    return Atom(position=np.zeros(3), mass_amu=35.0)


@pytest.fixture
def atom_h():
    return Atom(position=np.array([1.27455, 0.0, 0.0]), mass_amu=1.0)


def test_bond_length(atom_cl, atom_h):
    k_ab = 4.88e-8
    r_e_ab = 1.27455
    bond = Bond(atom_a=atom_cl, atom_b=atom_h, r_e_ab=r_e_ab, k_ab=k_ab)
    expected = np.array([r_e_ab, 0.0, 0.0])
    actual = bond.r_ab
    assert np.array_equal(expected, actual)
