using Test
using PyCall
using DFTK
bulk = pyimport("ase.build").bulk
include("calculator.jl")


py"""
import ase
cell = [[3.21, 0.0, 0.0], [-1.605, 2.7799415461480477, 0.0], [0.0, 0.0, 5.21304]]
positions = [[0.0, 0.0, 0.0], [0.0, 1.85329, 2.60652]]
atoms = ase.Atoms(symbols='Mg2', pbc=True, cell=cell, positions=positions)
"""

@testset "Basic LDA basis construction" begin
    calc = asedftk.DFTK(;xc="LDA", kpts=[1, 2, 3], ecut=300)
    calc.atoms = py"atoms"
    basis = calc.get_dftk_basis()

    lattice = basis.model.lattice
    @test lattice[:, 1] ≈ DFTK.units.Ǎ * [3.21, 0.0, 0.0]
    @test lattice[:, 2] ≈ DFTK.units.Ǎ * [-1.605, 2.7799415461480477, 0.0]
    @test lattice[:, 3] ≈ DFTK.units.Ǎ * [0.0, 0.0, 5.21304]

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

    xcterm = [tt for tt in basis.model.term_types if tt isa Xc][1]
    @test length(xcterm.functionals) == 1
    @test xcterm.functionals[1].identifier == :lda_xc_teter93

    @test length(basis.kpoints) == 4
    @test basis.Ecut ≈ 300 * DFTK.units.eV atol=1e-8
end

@testset "PBE semicore basis construction" begin
    calc = asedftk.DFTK(;xc="PBE", kpts=0.3, smearing=("Gaussian", 10),
                        pps="hgh.k")  # semicore psps
    calc.atoms = py"atoms"
    basis = calc.get_dftk_basis()

    atoms = basis.model.atoms
    @test length(atoms) == 1
    @test length(atoms) == 1
    @test all(at isa ElementPsp for (at, positions) in atoms)
    @test atoms[1][1].psp.identifier == "hgh/pbe/mg-q10.hgh"
    @test atoms[1][1].psp.Zion == 10

    @test basis.model.temperature ≈ 10 * DFTK.units.eV atol=1e-8
    @test basis.model.smearing isa DFTK.Smearing.Gaussian
    @test basis.model.spin_polarization == :none

    xcterm = [tt for tt in basis.model.term_types if tt isa Xc][1]
    @test length(xcterm.functionals) == 2
    @test xcterm.functionals[1].identifier == :gga_x_pbe
    @test xcterm.functionals[2].identifier == :gga_c_pbe

    kgrid = DFTK.kgrid_size_from_minimal_spacing(basis.model.lattice,
                                                 0.3 / DFTK.units.Ǎ)
    kcoords, ksymops = bzmesh_ir_wedge(kgrid, basis.model.symops)
    @test kcoords ≈ [kpt.coordinate for kpt in basis.kpoints]
end

@testset "Custom functional construction" begin
    kpts = [(0, 0, 0), (0.5, 0.7, 0.3)]
    calc = asedftk.DFTK(;xc="PBE", kpts=kpts, smearing=("Fermi-Dirac", 5))
    calc.atoms = py"atoms"
    basis = calc.get_dftk_basis()

    @test basis.model.temperature ≈ 5 * DFTK.units.eV atol=1e-8
    @test basis.model.smearing isa DFTK.Smearing.FermiDirac
    @test basis.model.spin_polarization == :none

    xcterm = [tt for tt in basis.model.term_types if tt isa Xc][1]
    @test length(xcterm.functionals) == 2
    @test xcterm.functionals[1].identifier == :gga_x_pbe
    @test xcterm.functionals[2].identifier == :gga_c_pbe

    @test length(basis.kpoints) == 2
    @test Float64.(basis.kpoints[1].coordinate) ≈ [kpts[1]...]
    @test Float64.(basis.kpoints[2].coordinate) ≈ [kpts[2]...]
end

@testset "Test kpoint options" begin
    kpointoptions = [
        (ase=0.5,                length=18, kpt1=[0.0, 0.0, 0.0], kpt2=[ 1/5, 0.0, 0.0]),
        (ase=[2, 3, 4],          length=12, kpt1=[0.0, 0.0, 0.0], kpt2=[-0.5, 0.0, 0.0]),
        (ase=[2, 3, 4, "gamma"], length=9,  kpt1=[1/4, 1/6, 0.0], kpt2=[-1/4, 1/6, 0.0]),
        (ase=[(0, 0.2, 0), (0.5, 0.2, 0.3)], length=2,
         kpt1=[0, 0.2, 0], kpt2=[0.5, 0.2, 0.3]),
    ]

    for params in kpointoptions
        calc = asedftk.DFTK(;kpts=params.ase)
        calc.atoms = py"atoms"
        basis = calc.get_dftk_basis()
        @test length(basis.kpoints) == params.length
        @test Float64.(basis.kpoints[1].coordinate) ≈ params.kpt1
        @test Float64.(basis.kpoints[2].coordinate) ≈ params.kpt2
    end
end

@testset "Test smearing options" begin
    smearingoptions = [
        (ase=("Fermi-Dirac", 10), temperature=10 * DFTK.units.eV,
         smearing=DFTK.Smearing.FermiDirac),
        (ase=("Gaussian", 5), temperature=5 * DFTK.units.eV,
         smearing=DFTK.Smearing.Gaussian),
        (ase=("Methfessel-Paxton", 5, 0), temperature=5 * DFTK.units.eV,
         smearing=DFTK.Smearing.Gaussian),
        (ase=("Methfessel-Paxton", 5, 1), temperature=5 * DFTK.units.eV,
         smearing=DFTK.Smearing.MethfesselPaxton1),
        (ase=("Methfessel-Paxton", 5, 2), temperature=5 * DFTK.units.eV,
         smearing=DFTK.Smearing.MethfesselPaxton2),
    ]

    for params in smearingoptions
        calc = asedftk.DFTK(;smearing=params.ase)
        calc.atoms = py"atoms"
        model = calc.get_dftk_model()
        @test model.smearing isa params.smearing
        @test model.temperature ≈ params.temperature atol=1e-8
    end
end

@testset "Test mixing options" begin
    mixingoptions = [
        (ase=("SimpleMixing"), mixing=SimpleMixing, α=1.0),
        (ase=("SimpleMixing", Dict("α" => 0.4, )), mixing=SimpleMixing, α=0.4),
        (ase=("KerkerMixing", Dict("α" => 0.4, "kF" => 0.7)),
         mixing=KerkerMixing, α=0.4),
        (ase=("KerkerMixing"), mixing=KerkerMixing, α=0.8),
    ]

    for params in mixingoptions
        calc = asedftk.DFTK(;mixing=params.ase)
        calc.atoms = py"atoms"
        mixing = calc.get_dftk_mixing()
        @test mixing isa params.mixing
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
    FORCES_PBE = [[0.0 0.0 0.0]; [0.0 0.0 0.0]]

    silicon = bulk("Si")
    silicon.calc = asedftk.DFTK(;xc="PBE", kpts=(3, 3, 3), ecut=190, scftol=1e-4,
                                nbands=12, label=label)
    @test silicon.get_potential_energy() ≈ ENERGY_PBE atol=1e-4 rtol=1e-4
    @test silicon.get_forces() ≈ FORCES_PBE atol=1e-2
    @test length(silicon.calc.scfres.eigenvalues[1]) ≥ 12

    orig_ene = silicon.get_potential_energy()
    orig_forces = silicon.get_forces()

    # Read resultsfile again:
    @test isfile(label * ".json")
    silicon.calc = asedftk.DFTK(restart=label)

    # Check we got the old parameters and results back:
    @test silicon.calc.parameters["xc"] == "PBE"
    @test silicon.calc.parameters["kpts"] == [3, 3, 3]
    @test silicon.calc.parameters["ecut"] == 190
    @test silicon.calc.results["energy"] ≈ orig_ene
    @test silicon.calc.results["forces"] ≈ orig_forces

    # Test some invariances
    @test silicon == asedftk.DFTK.read_atoms(label)
    @test silicon == asedftk.DFTK(restart=label, label=label).get_atoms()
end
