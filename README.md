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
1. Install Julia e.g. by [downloading the binary](https://julialang.org/downloads).
   The use of at least **Julia 1.4** is required.
   It is highly recommended you install Julia *before* installing asedftk.
2. Install asedftk from [PyPi](https://pypi.org/project/asedftk):
   ```
   pip install asedftk
   ```
   This automatically installs the [PyJulia](https://pypi.org/project/julia/) package,
   which allows Julia and Python codes to interoperate with each other.
3. Install the Julia dependencies of asedftk:
   ```
   python3 -c "import asedftk; asedftk.install()"
   ```
4. That's it, you're all set. But **please note**:
   Due to [some limitations](https://pyjulia.readthedocs.io/en/stable/troubleshooting.html#your-python-interpreter-is-statically-linked-to-libpython)
   in some Linux distros like Debian or Ubuntu
   you might need to run your Python scripts
   with the `python-jl` wrapper if you want to use asedftk in them.
   I.e. if you have written a calculation script `script.py` you
   might need to start it as `python-jl script.py`
   in order to be able to use asedftk.

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
