import asedftk

import pytest

from pytest import approx

from ase.io import read, write
from ase.build import bulk, molecule

#
# Adapted from the ASE test suite code in
# https://gitlab.com/ase/ase/-/blob/master/ase/test/calculator/test_al.py
#
@pytest.mark.skip("Collinear spin not yet in DFTK")
def test_aluminium_oxide():
    ENERGY_LDA = 0.0
    ENERGY_PBE = 0.0

    kpts = [2, 2, 2]
    params = {"scftol": 0.0001, "ecut": 200}

    calc = asedftk.DFTK(label="dftk", xc="LDA", kpts=kpts, **params)
    al = bulk("AlO", crystalstructure='rocksalt', a=4.5)
    al.calc = calc
    elda = al.get_potential_energy()
    assert elda == approx(ENERGY_LDA, rel=5e-3)

    calc.set(xc='PBE', kpts=kpts)
    epbe = al.get_potential_energy()
    assert not calc.calculation_required(al, ['energy'])
    assert epbe == approx(ENERGY_PBE, rel=5e-3)

    al = calc.get_atoms()
    elda2 = al.get_potential_energy()
    assert elda2 == approx(ENERGY_LDA, rel=5e-3)

    label = "dir/dftk-2"
    calc = asedftk.DFTK(label=label, atoms=al, xc="LDA", kpts=kpts, **params)
    elda3 = al.get_potential_energy()
    assert elda == elda3
    assert elda3 == approx(ENERGY_LDA, rel=5e-3)


#
# Adapted from the ASE test suite code in
# https://gitlab.com/ase/ase/-/blob/master/ase/test/calculator/test_h2.py
#
def test_hydrogen():
    ENERGY_H1_LDA = 0.0
    ENERGY_H1_PBE = 0.0
    ENERGY_H2_LDA = -29.56773881096904
    ENERGY_H2_PBE = -30.320772272053382
    BINIDNG_LDA = 2 * ENERGY_H1_LDA - ENERGY_H2_LDA
    BINIDNG_PBE = 2 * ENERGY_H1_PBE - ENERGY_H2_PBE

    kpts = [2, 2, 2]
    params = {"scftol": 0.0001, "ecut": 200}

    # H2 molecule
    calc = asedftk.DFTK(label="dftk", xc="LDA", kpts=kpts, **params)
    h2 = molecule('H2', calculator=calc, pbc=True)
    h2.center(vacuum=2.0)
    e2 = h2.get_potential_energy()
    calc.set(xc='PBE')
    e2pbe = h2.get_potential_energy()
    assert e2 == approx(ENERGY_H2_LDA, rel=5e-3)
    assert e2pbe == approx(ENERGY_H2_PBE, rel=5e-3)

    # TODO Remove this once we are there
    pytest.skip("DFTK does not yet have collinear spin.")

    # H atom
    h1 = h2.copy()
    del h1[1]
    h1.set_initial_magnetic_moments([1])
    h1.calc = calc
    e1pbe = h1.get_potential_energy()
    calc.set(xc='LDA')
    e1 = h1.get_potential_energy()
    assert e1 == approx(ENERGY_H1_LDA, rel=5e-3)
    assert e1pbe == approx(ENERGY_H1_PBE, rel=5e-3)

    try:
        m1 = h1.get_magnetic_moment()
    except NotImplementedError:
        pass
    else:
        assert m1 == 0.0

    # Check binding energy
    assert 2 * e1 - e2 == approx(BINIDNG_LDA, rel=5e-3)
    assert 2 * e1pbe - e2pbe == approx(BINIDNG_PBE, rel=5e-3)

    # Check some restart workflow
    calc = asedftk.DFTK(restart="dftk")
    assert calc.atoms == h1
    assert not calc.calculation_required(h1, ['energy'])
    h1 = calc.get_atoms()
    assert h1.get_potential_energy() == e1

    calc = asedftk.DFTK(label="dir/dftk-h1", atoms=h1, xc='LDA', kpts=kpts, **params)
    assert h1.get_potential_energy() == e1


#
# Adapted from the ASE test suite code in
# https://gitlab.com/ase/ase/-/blob/master/ase/test/calculator/test_traj.py
#
def test_trajectories():
    ENERGY = -30.38758246199375
    FORCES = [[0.0, 0.0,  1.88397764],
              [0.0, 0.0, -1.88397764]]

    h2 = molecule('H2', pbc=True)
    h2.center(vacuum=2.0)
    h2.calc = asedftk.DFTK()
    e = h2.get_potential_energy()
    assert not h2.calc.calculation_required(h2, ['energy'])
    assert e == approx(ENERGY, rel=5e-5, abs=1e-6)
    f = h2.get_forces()
    assert not h2.calc.calculation_required(h2, ['energy', 'forces'])
    for i in range(len(FORCES)):
        assert f[i] == approx(FORCES[i], rel=1e-2, abs=1e-3)

    # Writing and reading trajectories
    write('h2.traj', h2)
    h2 = read('h2.traj')
    assert abs(e - h2.get_potential_energy()) < 1e-12
    assert abs(f - h2.get_forces()).max() < 1e-12
