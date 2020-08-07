using LinearAlgebra
using FFTW
using PyCall

FFTW.set_num_threads(2)
BLAS.set_num_threads(2)

DFTKcalc = pyimport("asedftk").DFTK
ase_build = pyimport("ase.build")
BFGS = pyimport("ase.optimize").BFGS

# H2 with bond length 1 Ångström
h2 = ase_build.molecule("H2", pbc=true, vacuum=10)
h2.set_positions([[10 10 11.0]; [10 10 10]])
h2.calc = DFTKcalc(verbose=true, scftol=1e-6, xc="PBE", kpts=[1, 1, 1], ecut=50)

dyn = BFGS(h2, trajectory="H2.traj")
dyn.run(fmax=0.05)
