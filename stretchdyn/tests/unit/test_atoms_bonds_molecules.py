import pytest  # type: ignore
import numpy as np  # type: ignore
from stretchdyn.atoms_bonds import Atom, Bond, Molecule


@pytest.fixture
def hcl_equilibrium():
    r_e_ab = 1.27455
    k_ab = 4.88e-8
    cl = Atom(mass_amu=35, symbol="Cl", pos=np.zeros(3), vel=np.zeros(3))
    h = Atom(mass_amu=1, symbol="H", pos=np.array([r_e_ab, 0.0, 0.0]), vel=np.zeros(3))
    bond_cl_h = Bond(atom_a=cl, atom_b=h, r_e_ab=r_e_ab, k_ab=k_ab)
    bond_h_cl = Bond(atom_a=h, atom_b=cl, r_e_ab=r_e_ab, k_ab=k_ab)
    cl.bonds["h"] = bond_h_cl
    h.bonds["cl"] = bond_cl_h
    molecule = Molecule(atoms={"h": h, "cl": cl})
    return molecule


def test_hcl_r_ab(hcl_equilibrium):
    r_e_ab = 1.27455
    h = hcl_equilibrium.atoms["h"]
    cl = hcl_equilibrium.atoms["cl"]
    assert np.array_equal(h.bonds["cl"].r_ab, np.array([-r_e_ab, 0.0, 0.0]))
    assert np.array_equal(cl.bonds["h"].r_ab, np.array([r_e_ab, 0.0, 0.0]))


def test_hcl_unit_h_cl(hcl_equilibrium):
    expected = np.array([-1.0, 0.0, 0.0])
    h_cl_bond = hcl_equilibrium.atoms["h"].bonds["cl"]
    actual = h_cl_bond.unit
    assert np.array_equal(actual, expected)


def test_hcl_unit_cl_h(hcl_equilibrium):
    expected = np.array([1.0, 0.0, 0.0])
    cl_h_bond = hcl_equilibrium.atoms["cl"].bonds["h"]
    actual = cl_h_bond.unit
    assert np.array_equal(actual, expected)
