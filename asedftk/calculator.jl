using PyCall
import DFTK: PlaneWaveBasis, Model, self_consistent_field, load_psp
import DFTK: load_lattice, load_atoms, ElementPsp, model_DFT, forces
ase_units = pyimport("ase.units")
calculator = pyimport("ase.calculators.calculator")


#    self.label         contains place where results are written to
#    self.parameters    calculational parameters
#    self.results       result
@pydef mutable struct DFTK <: calculator.Calculator
    implemented_properties = ["energy", "free_energy", "forces"]

    function __init__(self; atoms=nothing, label="dftk", kwargs...)
        calculator.Calculator.__init__(self; label=label, kwargs...)

        # TODO Setup default_parameters and parameters in general
        self.default_parameters = Dict{String,Any}()
        self.scfres = nothing  # Results of the most recent SCF run
    end

    function reset(self)
        # Reset internal state by purging all intermediates
        calculator.Calculator.reset(self)
        self.scfres = nothing
    end

    function read(self, label)
        calculator.Calculator.read(label)
        # TODO Read atoms, parameters, results (and scfres) from json/hdf5 file
        println("Would read from json file " * label * ".json.")
    end

    function set(self; kwargs...)
        changed_parameters = calculator.Calculator.set(self; kwargs...)
        if !isempty(changed_parameters)
            # Currently reset on all parameter changes
            # TODO Improve
            self.reset()
        end
    end

    function check_state(self, atoms)
        system_changes = calculator.Calculator.check_state(self, atoms)
        if "pbc" in system_changes
            # Ignore boundary conditions (we always use PBCs)
            deleteat!(system_changes, findfirst(isequal("pbc"), system_changes))
        end
        system_changes
    end

    function get_dftk_model(self)
        # TODO Get from parameters
        functional = "lda"
        psp_family = "hgh"
        psp_core = "fullcore"

        # TODO Temperature, spin polarisation, smearing

        psploader(symbol) = load_psp(symbol, functional=functional,
                                     family=psp_family, core=Symbol(psp_core))


        lattice = load_lattice(self.atoms)
        atoms = [ElementPsp(element.symbol, psp=psploader(element.symbol)) => positions
                 for (element, positions) in load_atoms(self.atoms)]

        if functional == "lda"
            functionals = [:lda_xc_teter93]
        else
            functionals = [:gga_x_pbe, :gga_c_pbe]
        end
        model_DFT(lattice, atoms, functionals)
    end

    function get_dftk_basis(self; model=self.get_dftk_model())
        # TODO Get from parameters
        Ecut = 15
        kgrid = [1, 1, 1]
        PlaneWaveBasis(model, Ecut, kgrid=kgrid)
    end

    function calculate(self, atoms=nothing, properties=["energy"],
                       system_changes=calculator.all_changes)
        calculator.Calculator.calculate(self, atoms=atoms, properties=properties,
                                        system_changes=system_changes)
        results = self.results

        if isnothing(self.atoms)
            pyraise(calculator.CalculatorSetupError("An Atoms object must be provided to " *
                                                    "perform a calculation"))
        end

        # TODO Get from parameters
        tol = 1e-6
        scfres = self_consistent_field(self.get_dftk_basis(), tol=tol)

        results["energy"] = sum(values(scfres.energies)) * ase_units.Hartree

        # TODO Unit conversion!
        if "forces" in properties
            results["forces"] = forces(scfres)
        end

        # TODO Write results to disk
        self.results = results
    end
end
