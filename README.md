# DFTK-based calculator for ASE

| **Documentation**                                       | **Build Status**              |  **Installation**         |
|:------------------------------------------------------- |:----------------------------- |:------------------------- |
| [![][docs-img]][docs-url] [![][gitter-img]][gitter-url] | [![][travis-img]][travis-url] | [![][pypi-img]][pypi-url] |

[docs-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-url]: https://github.com/mfherbst/asedftk/blob/master/docs/asedftk.md

[gitter-img]: https://badges.gitter.im/DFTK-jl/community.svg
[gitter-url]: https://gitter.im/DFTK-jl/community

[travis-img]: https://api.travis-ci.com/mfherbst/asedftk.svg?branch=master
[travis-url]: https://travis-ci.com/mfherbst/asedftk

[pypi-img]: https://img.shields.io/pypi/v/asedftk
[pypi-url]: https://pypi.org/project/asedftk

Small wrapper around the
[**density-functional toolkit (DFTK)**](https://dftk.org)
to provide a
[calculator interface](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html)
compatible with [ASE](https://wiki.fysik.dtu.dk/ase/index.html),
the atomistic simulation environment.

DFTK is a small library of Julia algorithms
for experimentation with plane-wave-based
density-functional theory (PWDFT) methods.
Its performance is on the same order of magnitude as established packages
in the field. See [dftk.org](https://dftk.org) and
the [DFTK documentation](https://juliamolsim.github.io/DFTK.jl/dev/) for more details.

## Installation
See the [installation instructions](https://github.com/mfherbst/asedftk/blob/master/docs/installation.md).

## Basic usage
`asedftk.DFTK` is basically a class wrapping around DFTK and making it an
[ASE calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html).
Just use it like any other calculator class. For example:
```python
from asedftk import DFTK
from ase.build import bulk

atoms = bulk("Si")
atoms.calc = DFTK()
print(atoms.get_potential_energy())
```

More details can be found in the [asedftk documentation](https://github.com/mfherbst/asedftk/blob/master/docs/asedftk.md).
