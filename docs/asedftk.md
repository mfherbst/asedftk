# asedftk documentation

asedftk is a wrapper around the
[**density-functional toolkit (DFTK)**](https://dftk.org),
a Julia code for plane-wave based density-functional theory (PWDFT).
asedftk provides an [ASE](https://wiki.fysik.dtu.dk/ase/index.html)-compatible
calculator class,
which allows to directly employ DFTK inside ASE workflows.

## Installation
1. Install Julia e.g. by [downloading the binary](https://julialang.org/downloads).
   The use of at least Julia **1.5** is required.
1. Install asedftk from [PyPi](https://pypi.org/project/asedftk), for example
   using pip:
   ```
   pip install asedftk
   ```
   Note: You can also use `pip` if you are managing your python packages
   with anaconda (`conda install pip` followed by above command), but it
   is not recommended to do this in the root anaconda environment.

This will automatically install DFTK inside a separate Julia environment
private to asedftk. Typically this DFTK version will be updated as needed.
In case you need to do the update manually at some point, simply run
`asedftk.update()`.

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
	- `"SimpleMixing(α=0.7)"`: Simple mixing with damping `0.7`
	- `"KerkerMixing(α=0.7, kTF=1.0)"`: Kerker mixing with damping `α = 0.7` and
	  screening parameter `kTF = 1.0`. Applies the operator kernel
	  `α * G^2 / (kTF^2 + G^2)` for a wave vector `G` in frequency space.
	- `"HybridMixing(α=0.7, kTF=1.0, εr=10)"`: LDOS mixing as described
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
- **xc**: Exchange-correlation functional, default is `LDA`, options are `LDA` or `PBE`.
- **n_mpi**: Number of MPI threads to employ when running DFTK.
  Set this to a value between `1` and the number of irreducible k-Points.
  See the [DFTK parallelization documentation](https://juliamolsim.github.io/DFTK.jl/stable/guide/parallelization/) for details.
- **n_threads**: Number of Julia threads to employ when running DFTK.
  Typically the best setting is to use only one thread (the default).
  See the [DFTK parallelization documentation](https://juliamolsim.github.io/DFTK.jl/stable/guide/parallelization/) for details.


## Custom Julia environment for asedftk
If you want to employ a custom Julia environment in combination with asedftk,
for example to supply a modified DFTK,
just point the environment variable `ASEDFTK_DFTK_ENVIRONMENT` to the folder with your
custom environment.
When doing this, ensure that the packages specified in asedftk's
[Project.toml](../asedftk/dftk_environment/Project.toml)
are available in your custom environment.
