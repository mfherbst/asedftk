# DFTK-based calculator for ASE
[![Documentation](https://img.shields.io/badge/doc-latest-blue.svg)](https://github.com/mfherbst/asedftk/blob/master/docs/asedftk.md)
[![pypi](https://img.shields.io/pypi/v/asedftk)](https://pypi.org/project/asedftk)
[![gitter](https://badges.gitter.im/DFTK-jl/community.svg)](https://gitter.im/DFTK-jl/community)
[![license](https://img.shields.io/github/license/mfherbst/asedftk.svg?maxAge=2592000)](https://github.com/mfherbst/asedftk/blob/master/LICENSE)
[![Build Status](https://api.travis-ci.com/mfherbst/asedftk.svg?branch=master)](https://travis-ci.com/mfherbst/asedftk)

Small wrapper around the
[**density-functional toolkit (DFTK)**](https://dftk.org)
to provide a
[calculator interface](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html)
compatible with [ASE](https://wiki.fysik.dtu.dk/ase/index.html),
the atomistic simulation environment.

DFTK is a small library of Julia algorithms
for experimentation with plane-wave-based
density-functional theory (PWDFT) methods.
While this code is not yet fully optimised and fine-tuned,
performance is on the same order of magnitude as established packages
in the field. See [dftk.org](https://dftk.org) for more details.

## Installation
1. Install Julia e.g. by [downloading the binary](https://julialang.org/downloads).
   The use of at least **Julia 1.4** is required.
   It is highly recommended you do this **before** installing asedftk.
2. Install asedftk from [PyPi](https://pypi.org/project/asedftk):
   ```
   pip install asedftk
   ```
   This will also install the [PyJulia](https://pypi.org/project/julia/) package,
   which allows Julia and Python to interoperate.
3. Install Julia dependencies:
   ```
   python-jl -c "import asedftk; asedftk.install()"
   ```
   **Note:** The use of `python-jl` instead of a plain `python`
   is on purpose
   to [work around some limitations](https://pyjulia.readthedocs.io/en/stable/troubleshooting.html#your-python-interpreter-is-statically-linked-to-libpython)
   present in some Linux distros like Debian or Ubuntu.

## Basic usage
`asedftk.DFTK` is basically a class wrapping around DFTK and making it an
[ASE calculator](https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html).
Just use it like any other. For example:
```python
from asedftk import DFTK
from ase.build import bulk

atoms = bulk("Si")
atoms.calc = DFTK()
print(atoms.get_potential_energy())
```
Keep in mind that in order to use this you might again need
to call this script with the `python-jl` wrapper
instead of just using `python` as usual.

More details can be found in the [asedftk documentation](https://github.com/mfherbst/asedftk/blob/master/docs/asedftk.md)
