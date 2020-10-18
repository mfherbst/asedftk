# asedftk documentation

asedftk is a wrapper around the
[**density-functional toolkit (DFTK)**](https://dftk.org),
a Julia code for plane-wave based density-functional theory (PWDFT).
asedftk provides an [ASE](https://wiki.fysik.dtu.dk/ase/index.html)-compatible
calculator class,
which allows to directly employ DFTK inside ASE workflows.

## Installation
See the [installation instructions](installation.md).

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
magnesium.calc = DFTK(xc="PBE", smearing=("Gaussian", 0.027), kpts=(5, 5, 5))
print("Magnesium energy: ", magnesium.get_potential_energy())
print("Magnesium forces: ", magnesium.get_forces())
```

### Hydrogen geometry optimisation
Using the ASE optimiser from python:
```python
# Python
from asedftk import DFTK
import ase.build
from ase.optimize import BFGS

h2 = ase.build.molecule("H2", pbc=True, vacuum=10)
h2.set_positions([[10, 10, 11.0], [10, 10, 10]])
h2.calc = DFTK(scftol=1e-6, xc="PBE", kpts=[1, 1, 1], ecut=50)

dyn = BFGS(h2, trajectory="H2.traj")
dyn.run(fmax=0.05)

print("H-H distance: ", h2.get_distance(0, 1))
```
Alternatively using a Julia script (and `PyCall`) to drive the optimisation:
```julia
# Julia
using PyCall
DFTKcalc = pyimport("asedftk").DFTK
ase_build = pyimport("ase.build")
BFGS = pyimport("ase.optimize").BFGS

h2 = ase_build.molecule("H2", pbc=true, vacuum=10)
h2.set_positions([[10 10 11.0]; [10 10 10]])
h2.calc = DFTKcalc(scftol=1e-6, xc="PBE", kpts=[1, 1, 1], ecut=50)

dyn = BFGS(h2, trajectory="H2.traj")
dyn.run(fmax=0.05)

println("H-H distance: ", h2.get_distance(0, 1))
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
	- `(n1,n2,n3)`: (Unshifted) Monkhorst-Pack grid
	- `(n1,n2,n3,"gamma")`: Shifted MP grid to contain the gamma point.
	- `[(k11,k12,k13),(k21,k22,k23),...]`: Explicit k-Point list in units of the reciprocal lattice vectors
    - `3.5` (or any float): k-point density as in `3.5` kpoints per Ǎngström.
- **mixing**: Mixing scheme used during SCF iterations. Examples for valid options:
	- `("SimpleMixing", dict(α=0.7))`: Simple mixing with damping `0.7`
	- `("KerkerMixing", dict(α=0.7, kTF=1.0))`: Kerker mixing with damping `α = 0.7` and
	  screening parameter `kTF = 1.0`. Applies the operator kernel
	  `α * G^2 / (kTF^2 + G^2)` for a wave vector `G` in frequency space.
	- `("HybridMixing", dict(α=0.7, kTF=1.0, εr=10)`: LDOS mixing as described
	  in [arXiv 2009.01665](https://arxiv.org/abs/2009.01665).
	  If `εr=1` a plain LDOS mixing is used, for other values a model dielectric
	  term is added as well. `kTF` is meaningless unless the model dielectric term
	  is used.
- **nbands**: Number of bands to compute
- **pps**: Pseudopotential family. Currently the only choices are `"hgh"`
  (Goedecker-type pseudos) and `"hgh.k"` (semi-core version of `"hgh"`).
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
  For details on the parameters, which influence threading in DFTK,
  see [the DFTK documentation](https://docs.dftk.org/dev/guide/parallelisation/).
  To set BLAS and FFTW threads from python use:
  ```
  from julia.LinearAlgebra import BLAS
  from julia import FFTW

  FFTW.set_num_threads(N)
  BLAS.set_num_threads(N)
  ```
- On Debian and Ubuntu the use of the `python-jl` wrapper script
  is unfortunately necessary due to the `python` executable being
  statically linked to the libpython.
  See the [PyJulia documentation](https://pyjulia.readthedocs.io/en/stable/troubleshooting.html#your-python-interpreter-is-statically-linked-to-libpython)
  for details.
- There are a few more [smaller limitations](https://pyjulia.readthedocs.io/en/stable/limitations.html)
  in the interoperation of Python and Julia you probably should give a read to be aware.
