from numpy.testing import assert_allclose

import asedftk

from ase.build import bulk

ENERGY_PBE = -213.12688268374683  # eV
FORCES_PBE = [[25.34636840001346,  9.7064846689645, -3.361631681010807e-7],
              [-25.34636754014557, -9.7064849446649, -3.689230622334620e-7]]


def test_silicon():
    silicon = bulk("Si")
    silicon.calc = asedftk.DFTK(xc="PBE", kpts=(3, 3, 3), ecut=190, scftol=1e-4)
    assert_allclose(silicon.get_potential_energy(), ENERGY_PBE, atol=1e-4)
    assert_allclose(silicon.get_forces(), FORCES_PBE, atol=5e-3)
