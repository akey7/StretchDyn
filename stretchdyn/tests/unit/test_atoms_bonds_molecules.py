import pytest  # type: ignore
import numpy as np  # type: ignore
from stretchdyn.atoms_bonds import Atom, Bond, Molecule


@pytest.fixture
def hcl_equilibrium():
    r_e_ab = 1.27455
    k_ab = 4.88e-8
    cl = Atom(
        mass_amu=35, symbol="Cl", pos=np.zeros(3), vel=np.zeros(3), prev_accel=np.zeros(3)
    )
    h = Atom(
        mass_amu=1,
        symbol="H",
        pos=np.array([r_e_ab, 0.0, 0.0]),
        vel=np.zeros(3),
        prev_accel=np.zeros(3),
    )
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


def test_hcl_stretch_derivative_h_cl(hcl_equilibrium):
    expected = 0.0
    h_cl_bond = hcl_equilibrium.atoms["h"].bonds["cl"]
    actual = h_cl_bond.stretch_derivative
    assert actual == expected


def test_hcl_stretch_derivative_cl_h(hcl_equilibrium):
    expected = 0.0
    cl_h_bond = hcl_equilibrium.atoms["cl"].bonds["h"]
    actual = cl_h_bond.stretch_derivative
    assert actual == expected


def test_hcl_stretch_force_h_cl(hcl_equilibrium):
    expected = np.array([0.0, 0.0, 0.0])
    h_cl_bond = hcl_equilibrium.atoms["h"].bonds["cl"]
    actual = h_cl_bond.stretch_force
    assert np.array_equal(actual, expected)


def test_hcl_stretch_force_cl_h(hcl_equilibrium):
    expected = np.array([0.0, 0.0, 0.0])
    cl_h_bond = hcl_equilibrium.atoms["cl"].bonds["h"]
    actual = cl_h_bond.stretch_force
    assert np.array_equal(actual, expected)


def test_hcl_net_stretch_force_h(hcl_equilibrium):
    expected = np.array([0.0, 0.0, 0.0])
    actual = hcl_equilibrium.atoms["h"].net_stretch_force
    assert np.array_equal(actual, expected)


def test_hcl_net_stretch_force_cl(hcl_equilibrium):
    expected = np.array([0.0, 0.0, 0.0])
    actual = hcl_equilibrium.atoms["cl"].net_stretch_force
    assert np.array_equal(actual, expected)


def test_hcl_net_acceleration_h(hcl_equilibrium):
    expected = np.array([0.0, 0.0, 0.0])
    actual = hcl_equilibrium.atoms["h"].stretch_accel
    assert np.array_equal(actual, expected)


def test_hcl_net_acceleration_cl(hcl_equilibrium):
    expected = np.array([0.0, 0.0, 0.0])
    actual = hcl_equilibrium.atoms["cl"].stretch_accel
    assert np.array_equal(actual, expected)


def test_hcl_update_pos_h(hcl_equilibrium):
    r_e_ab = 1.27455
    expected = np.array([r_e_ab, 0.0, 0.0])
    h = hcl_equilibrium.atoms["h"]
    h.update_pos_vel()
    actual = h.pos
    assert np.array_equal(actual, expected)


def test_hcl_update_pos_cl(hcl_equilibrium):
    expected = np.array([0.0, 0.0, 0.0])
    cl = hcl_equilibrium.atoms["cl"]
    cl.update_pos_vel()
    actual = cl.pos
    assert np.array_equal(actual, expected)


def test_hcl_update_vel_h(hcl_equilibrium):
    expected = np.array([0, 0.0, 0.0])
    h = hcl_equilibrium.atoms["h"]
    h.update_pos_vel()
    actual = h.vel
    assert np.array_equal(actual, expected)


def test_hcl_update_vel_cl(hcl_equilibrium):
    expected = np.array([0.0, 0.0, 0.0])
    cl = hcl_equilibrium.atoms["cl"]
    cl.update_pos_vel()
    actual = cl.vel
    assert np.array_equal(actual, expected)
