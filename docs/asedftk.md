# asedftk documentation

asedftk is a wrapper around the
[**density-functional toolkit (DFTK)**](https://dftk.org),
a Julia code for plane-wave based density-functional theory (PWDFT).
asedftk provides an [ASE](https://wiki.fysik.dtu.dk/ase/index.html)-compatible
calculator class,
which allows to directly employ DFTK inside ASE workflows.

## Installation
See the [README](../README.md#installation).

## Some examples
The best way to explain asedftk is by some examples:

### Silicon energy and forces
```python
from asedftk import DFTK
import ase.build

silicon = ase.build.bulk("Si")
silicon.calc = DFTK()

print("Silicon energy: ", silicon.get_potential_energy())
print("Silicon forces: ", silicon.get_forces())
```

### Magnesium energy and forces
```python
from asedftk import DFTK
import ase.build

magnesium = ase.build.bulk("Mg")
magnesium.calc = DFTK(xc="PBE", smearing=("Gaussian", 10), nbands=8, kpts=(5, 5, 5))
print("Magnesium energy: ", magnesium.get_potential_energy())
print("Magnesium forces: ", magnesium.get_forces())
```

## DFTK parameters
The `asedftk.DFTK` class supports a couple of parameters
to customise the calculation,
mostly following the
[ASE calculator interface proposal](https://wiki.fysik.dtu.dk/ase/development/proposals/calculators.html).
See also the `__init__` function
in [calculator.jl](https://github.com/mfherbst/asedftk/blob/master/asedftk/calculator.jl).

- **ecut**: Kinetic energy cutoff. 400eV by default.
- **functionals**: Explicit functional list.
  This overwrites the choice of **xc**, but **xc** still needs to be given to select
  appropriate pseudopotentials. Valid options are
  all [functionals from libxc](https://www.tddft.org/programs/libxc/functionals/),
  for example `"lda_x"`, `"GGA_X_B86"`, `"GGA_C_LYP"`.
- **kpts**: k-point grid to use. Valid options are:
	- `(n1,n2,n3)`: Monkhorst-Pack grid
	- `(n1,n2,n3,"gamma")`: Shifted MP grid to contain the gamma point.
	- `[(k11,k12,k13),(k21,k22,k23),...]`: Explicit k-Point list in units of the reciprocal lattice vectors
	- `3.5` (or any float): k-point density as in `3.5` kpoints per Ǎngström.
- **mixing**: Mixing scheme used during SCF iterations. Examples for valid options:
	- `("SimpleMixing", "0.7")`: Simple mixing with damping `0.7`
	- `("KerkerMixing", "0.7", "1.0")`: Kerker mixing with damping `α=0.7` and
	  screening parameter `G0=1.0`. Applies the operator kernel
	  `α * G^2 / (G0^2 + G^2)` for a wave vector `G` in frequency space.
- **nbands**: Number of bands to compute
- **pps**: Pseudopotential family. Currently the only choices are `hgh`
  (Goedecker-type pseudos) and `hgh.k` (semi-core version of `hgh`).
- **scftol**: Convergence tolerance of the SCF. Default `1e-5`.
- **smearing**: Smearing function and temperature, options:
	- `('Fermi-Dirac', width)`
	- `('Gaussian', width)`
	- `('Methfessel-Paxton', width, n)`
	where in each case `width` is the width in eV and `n` is the Methfessel-Paxton order.
- **verbose**: Make the SCF be more verbose and print some iteration information.
- **xc**: Exchange-correlation functional, default is `LDA`, options are `LDA` or `PBE`.

## Tips and tricks
Generally combining Python and Julia codes (like this package does) is seamless.
Still there are a few rough edges,
which are worth knowing when working with asedftk:

- **Threading:** DFTK makes use of Julia's threading framework.
  This means that by default it runs single-threaded.
  To benefit from a multi-core CPU you need to explicitly set
  the environment variable `JULIA_NUM_THREADS`
  to a value above 1 *before* using asedftk,
  i.e. before even starting the script.
- On Debian and Ubuntu the use of the `python-jl` wrapper script
  is unfortunately necessary due to the `python` executable being
  statically linked to the libpython.
  See the [PyJulia documentation](https://pyjulia.readthedocs.io/en/stable/troubleshooting.html#your-python-interpreter-is-statically-linked-to-libpython)
  for details.
- There are a few more [smaller limitations](https://pyjulia.readthedocs.io/en/stable/limitations.html)
  in the interoperation of Python and Julia you probably should give a read to be aware.
