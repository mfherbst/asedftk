from numpy.testing import assert_allclose

import asedftk

from ase.build import bulk

ENERGY_PBE = -854.8675691098085  # eV
FORCES_PBE = [[-6.316188006537352, -6.316035754751137, -14.510811801671302],
              [6.310114551209851, 6.309725717873144, 14.504810090181293],
              [-6.316527336116248, -6.31637644522554, -14.50856193248493],
              [6.311023337236977, 6.322789139759234, 14.51507332513539],
              [-6.3162593328581105, -6.316409108548032, -14.508542464053221],
              [6.3223038282193516, 6.310403017560323, 14.515101406186925],
              [-6.316312396627089, -6.316506421049234, -14.511626399173634],
              [6.322487450255041, 6.323020322615514, 14.505606155836926]]


def test_potential_energy():
    silicon = bulk("Si", cubic=True)
    silicon.calc = asedftk.DFTK(xc="PBE", kpts=(3, 3, 3), ecut=190, scftol=1e-4)
    assert_allclose(silicon.get_potential_energy(), ENERGY_PBE, atol=1e-4)
    assert_allclose(silicon.get_forces(), FORCES_PBE, atol=5e-3)
