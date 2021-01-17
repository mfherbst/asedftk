using Test
using JSON
using LinearAlgebra
using Unitful
using UnitfulAtomic
include("run_calculation.jl")

@testset "Basic LDA basis construction" begin
    state      = load_state("mg.json")
    parameters = Dict(state["parameters"]...,
                      "xc"   => "LDA",
                      "kpts" => [1, 2, 3],
                      "ecut" => 300)
    basis = get_dftk_basis(parameters, state["extra"])

    lattice = basis.model.lattice
    @test lattice[:, 1] ≈ austrip(1u"Å") * [3.21, 0.0, 0.0]
    @test lattice[:, 2] ≈ austrip(1u"Å") * [-1.605, 2.7799415461480477, 0.0]
    @test lattice[:, 3] ≈ austrip(1u"Å") * [0.0, 0.0, 5.21304]

    atoms = basis.model.atoms
    @test length(atoms) == 1
    @test all(at isa ElementPsp for (at, positions) in atoms)
    @test atoms[1][1].symbol == :Mg
    @test atoms[1][2][1] ≈ [0, 0, 0]
    @test atoms[1][2][2] ≈ [1/3, 2/3, 1/2] atol=1e-5
    @test atoms[1][1].psp.identifier == "hgh/lda/mg-q2.hgh"
    @test atoms[1][1].psp.Zion == 2

    @test basis.model.temperature == 0.0
    @test basis.model.smearing isa DFTK.Smearing.None
    @test basis.model.spin_polarization == :none

    xcterm = only(tt for tt in basis.model.term_types if tt isa Xc)
    @test only(xcterm.functionals) == :lda_xc_teter93

    @test length(basis.kpoints) == 4
    @test basis.Ecut ≈ 300 * austrip(1u"eV") atol=1e-8
end

@testset "PBE semicore basis construction" begin
    state      = load_state("mg.json")
    parameters = Dict(state["parameters"]...,
                      "xc"       => "PBE",
                      "kpts"     => 3.333,
                      "smearing" => ("Gaussian", 10),
                      "pps"      => "hgh.k")  # semicore psps
    basis = get_dftk_basis(parameters, state["extra"])

    atoms = basis.model.atoms
    @test length(atoms) == 1
    @test length(atoms) == 1
    @test all(at isa ElementPsp for (at, positions) in atoms)
    @test atoms[1][1].psp.identifier == "hgh/pbe/mg-q10.hgh"
    @test atoms[1][1].psp.Zion == 10

    @test basis.model.temperature ≈ 10 * austrip(1u"eV") atol=1e-8
    @test basis.model.smearing isa DFTK.Smearing.Gaussian
    @test basis.model.spin_polarization == :none

    xcterm = only(tt for tt in basis.model.term_types if tt isa Xc)
    @test length(xcterm.functionals) == 2
    @test xcterm.functionals[1] == :gga_x_pbe
    @test xcterm.functionals[2] == :gga_c_pbe

    kgrid = DFTK.kgrid_size_from_minimal_spacing(basis.model.lattice,
                                                 1 / 3.333austrip(1u"Å"))
    kcoords, ksymops = bzmesh_ir_wedge(kgrid, basis.model.symmetries)
    @test kcoords ≈ [kpt.coordinate for kpt in basis.kpoints]
end

@testset "Custom functional construction" begin
    kpts = [(0, 0, 0), (0.5, 0.7, 0.3)]
    state      = load_state("mg.json")
    parameters = Dict(state["parameters"]...,
                      "xc"       => "PBE",
                      "kpts"     => kpts,
                      "smearing" => ("Fermi-Dirac", 5))
    basis = get_dftk_basis(parameters, state["extra"])

    @test basis.model.temperature ≈ 5 * austrip(1u"eV") atol=1e-8
    @test basis.model.smearing isa DFTK.Smearing.FermiDirac
    @test basis.model.spin_polarization == :none

    xcterm = only(tt for tt in basis.model.term_types if tt isa Xc)
    @test length(xcterm.functionals) == 2
    @test xcterm.functionals[1] == :gga_x_pbe
    @test xcterm.functionals[2] == :gga_c_pbe

    @test length(basis.kpoints) == 2
    @test Float64.(basis.kpoints[1].coordinate) ≈ [kpts[1]...]
    @test Float64.(basis.kpoints[2].coordinate) ≈ [kpts[2]...]
end

@testset "Test kpoint options" begin
    kpointoptions = [
        (ase=2.0,                length=18, kpt1=[0.0, 0.0, 0.0], kpt2=[ 1/5, 0.0, 0.0]),
        (ase=[2, 3, 4],          length=12, kpt1=[0.0, 0.0, 0.0], kpt2=[-0.5, 0.0, 0.0]),
        (ase=[2, 3, 4, "gamma"], length=6,  kpt1=[1/4, 0.0, 1/8], kpt2=[ 1/4, 1/3, 1/8]),
        (ase=[(0, 0.2, 0), (0.5, 0.2, 0.3)], length=2,
         kpt1=[0, 0.2, 0], kpt2=[0.5, 0.2, 0.3]),
    ]

    for params in kpointoptions
        state      = load_state("mg.json")
        parameters = Dict(state["parameters"]..., "kpts" => params.ase)
        basis = get_dftk_basis(parameters, state["extra"])
        @test length(basis.kpoints) == params.length
        @test Float64.(basis.kpoints[1].coordinate) ≈ params.kpt1
        @test Float64.(basis.kpoints[2].coordinate) ≈ params.kpt2
    end
end

@testset "Test smearing options" begin
    smearingoptions = [
        (ase=("Fermi-Dirac", 10), temperature=10 * austrip(1u"eV"),
         smearing=DFTK.Smearing.FermiDirac),
        (ase=("Gaussian", 5), temperature=5 * austrip(1u"eV"),
         smearing=DFTK.Smearing.Gaussian),
        (ase=("Methfessel-Paxton", 5, 0), temperature=5 * austrip(1u"eV"),
         smearing=DFTK.Smearing.Gaussian),
        (ase=("Methfessel-Paxton", 5, 1), temperature=5 * austrip(1u"eV"),
         smearing=DFTK.Smearing.MethfesselPaxton1),
        (ase=("Methfessel-Paxton", 5, 2), temperature=5 * austrip(1u"eV"),
         smearing=DFTK.Smearing.MethfesselPaxton2),
    ]

    for params in smearingoptions
        state      = load_state("mg.json")
        parameters = Dict(state["parameters"]..., "smearing" => params.ase)
        model = get_dftk_model(parameters, state["extra"])
        @test model.smearing isa params.smearing
        @test model.temperature ≈ params.temperature atol=1e-8
    end
end

@testset "Test mixing options" begin
    mixingoptions = [
        (ase="SimpleMixing()",               mixing=SimpleMixing, α=0.8),
        (ase="SimpleMixing(α=0.4)",          mixing=SimpleMixing, α=0.4),
        (ase="KerkerMixing(α=0.4, kTF=0.7)", mixing=KerkerMixing, α=0.4),
        (ase="KerkerMixing()",               mixing=KerkerMixing, α=0.8),
        (ase="DielectricMixing()",           mixing=DielectricMixing, α=0.8),
        (ase="DielectricMixing(α=0.2, kTF=0.7, εr=8.0)",
                                             mixing=DielectricMixing, α=0.2),
        # HybridMixing gets translated to DFTK.χ0Mixing internally
        (ase="HybridMixing()",               mixing=χ0Mixing, α=0.8),
        (ase="HybridMixing(α=0.2, kTF=0.7, εr=8.0)",
                                             mixing=χ0Mixing, α=0.2),
    ]

    for params in mixingoptions
        state      = load_state("mg.json")
        parameters = Dict(state["parameters"]..., "mixing" => params.ase)
        mixing = get_dftk_mixing(parameters, state["extra"])
        @test (mixing isa params.mixing)
        @test mixing.α == params.α
    end
end

# TODO Missing checks for mandatory ASE parameters:
#      See https://wiki.fysik.dtu.dk/ase/development/calculators.html
#      and https://wiki.fysik.dtu.dk/ase/development/proposals/calculators.html
#      - charge

@testset "Silicon calculation" begin
    label = "silicon_dftk"
    ENERGY_PBE = -213.12688268374683  # eV
    FORCES_PBE = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]

    state = load_state("si.json")
    state["parameters"] = Dict(state["parameters"]...,
                               "xc" => "PBE",
                               "kpts" => (3, 3, 3),
                               "ecut" => 190,
                               "scftol" => 1e-4,
                               "nbands" => 12)
    save_state(label * ".json", state)

    rm(state["scfres"], force=true)
    state = run_calculation(["energies", "forces"], label * ".json")
    @test state["results"]["energy"] ≈ ENERGY_PBE atol=1e-4 rtol=1e-4
    @test state["results"]["forces"] ≈ FORCES_PBE atol=1e-2

    orig_ene    = state["results"]["energy"]
    orig_forces = state["results"]["forces"]
    orig_atoms  = state["atoms"]

    save_state(label * ".json", state)
    state = load_state(label * ".json")

    # Check we got the old parameters and results back:
    @test state["parameters"]["xc"]   == "PBE"
    @test state["parameters"]["kpts"] == [3, 3, 3]
    @test state["parameters"]["ecut"] == 190
    @test state["results"]["energy"] ≈ orig_ene
    @test state["results"]["forces"] ≈ orig_forces
    @test state["atoms"] == orig_atoms
end
