using DFTK
using JLD2
using JSON

# DFTK
setup_threading()

a = 10.26  # Silicon lattice constant in Bohr
lattice = a / 2 * [[0 1 1.];
                   [1 0 1.];
                   [1 1 0.]]
Si = ElementPsp(:Si, psp=load_psp("hgh/lda/Si-q4"))
atoms = [Si => [ones(3)/8, -ones(3)/8]]

model = model_LDA(lattice, atoms)
basis = PlaneWaveBasis(model, 10.0; kgrid=[1, 1, 1])

callback = DFTK.ScfSaveCheckpoints("checkpointfile.jld2") ∘ DFTK.ScfDefaultCallback()
scfres = self_consistent_field(basis, tol=1e-2, callback=callback, mixing=KerkerMixing())
compute_forces_cart(scfres)

# JSON
JSON.parse(JSON.json(Dict("a" => 1.0)))
