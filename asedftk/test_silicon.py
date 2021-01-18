import asedftk

from numpy.testing import assert_allclose

from ase.build import bulk

ENERGY_PBE = -213.12688268374683  # eV
FORCES_PBE = [[0, 0, 0.], [0, 0, 0.]]


def test_silicon():
    silicon = bulk("Si")
    silicon.calc = asedftk.DFTK(xc="PBE", kpts=(3, 3, 3), ecut=190, scftol=1e-6,
                                mixing="SimpleMixing(α=1.3)")
    assert_allclose(silicon.get_potential_energy(), ENERGY_PBE, atol=1e-5)
    assert_allclose(silicon.get_forces(), FORCES_PBE, atol=1e-2)
