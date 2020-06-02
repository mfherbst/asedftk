module asedftk
using PyCall
using LinearAlgebra
import JSON
import DFTK: PlaneWaveBasis, Model, self_consistent_field, load_psp
import DFTK: load_lattice, load_atoms, ElementPsp, model_DFT, forces
import DFTK: Smearing, kgrid_size_from_minimal_spacing, Vec3, Mat3
import DFTK: KerkerMixing, SimpleMixing, scf_default_callback

ase_units = pyimport("ase.units")
ase_io = pyimport("ase.io")
calculator = pyimport("ase.calculators.calculator")

# Raise an InputError exception
inputerror(s) = pyraise(calculator.InputError(s))

# Discussion of valid parameters:
#
# - xc, kpts, smearing, nbands as in
#   https://wiki.fysik.dtu.dk/ase/development/proposals/calculators.html
# - functional: If specified, defines the precise functionals used from libxc,
#   if not autodetermined from xc parameter.
# - pps: Pseudopotential family employed. Valid are "hgh" and "hgh.k" in which
#   case the semi-core versions of the HGH pseudopotentials is used.
# - scftol: SCF convergence tolerance.
# - ecut: Kinetic energy cutoff (in eV)
#

@pydef mutable struct DFTK <: calculator.Calculator
    implemented_properties = ["energy", "forces"]

    function __init__(self; atoms=nothing, label="dftk", kwargs...)
        self.default_parameters = Dict{String,Any}(
            "xc"          => "LDA",
            "kpts"        => 0.25,
            "smearing"    => nothing,
            "nbands"      => nothing,
            "charge"      => 0.0,
            "functionals" => nothing,
            "pps"         => "hgh",
            "scftol"      => 1e-5,
            "ecut"        => 400,
            "verbose"     => false,
            "mixing"      => nothing,
        )

        calculator.Calculator.__init__(self; label=label, kwargs...)
        self.scfres = nothing   # Results of the most recent SCF run
    end

    function reset(self)
        # Reset internal state by purging all intermediates
        calculator.Calculator.reset(self)
        self.scfres = nothing
    end

    function read(self, label)
        calculator.Calculator.read(self, label)

        saved_dict = open(JSON.parse, label * ".json", "r")
        self.parameters = calculator.Parameters(saved_dict["parameters"])
        self.results = saved_dict["results"]
        @pywith pyimport("io").StringIO(saved_dict["atoms"]) as f begin
            self.atoms = ase_io.read(f, format="json")
        end

        self.scfres = nothing  # Discard any cached results
    end

    function write(self, label=self.label)
        if isnothing(self.atoms)
            pyraise(calculator.CalculatorSetupError("An Atoms object must be present " *
                                                    "in the calculator to use write"))
        end

        # Convert atoms to json
        atoms_json = ""
        @pywith pyimport("io").StringIO() as f begin
            ase_io.write(f, self.atoms.copy(), format="json")
            atoms_json = f.getvalue()
        end

        # Transpose arrays inside results, such that they are written
        # in a predictable order to disk
        results = Dict{String, Any}()
        for (key, val) in pairs(self.results)
            if val isa AbstractArray
                if ndims(val) == 2
                    results[key] = collect(eachrow(val))
                    continue
                elseif ndims(val) > 2
                    error("Arrays of dimension larger 2 not supported")
                end
            end
            results[key] = val
        end
        open(label * ".json", "w") do fp
            save_dict = Dict(
                "parameters" => self.parameters,
                "results" => results,
                "atoms" => atoms_json,
            )
            JSON.print(fp, save_dict)
        end
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

    function get_spin_polarized(self)
        if isnothing(self.scfres)
            false
        else
            self.scfres.basis.model.spin_polarization in (:full, :collinear)
        end
    end

    function get_dftk_model(self)
        inputerror_param(param) = inputerror("Unknown value to $param: " *
                                             "$(self.parameters[param])")
        if isnothing(self.atoms)
            pyraise(calculator.CalculatorSetupError("An Atoms object must be provided to " *
                                                    "perform a calculation"))
        end

        # Parse psp and DFT functional parameters
        functionals = [:lda_xc_teter93]
        if lowercase(self.parameters["xc"]) == "lda"
            psp_functional = "lda"
            functionals = [:lda_xc_teter93]
        elseif lowercase(self.parameters["xc"]) == "pbe"
            psp_functional = "pbe"
            functionals = [:gga_x_pbe, :gga_c_pbe]
        else
            inputerror_param("xc")
        end

        if !isnothing(self.parameters["functionals"])
            funs = self.parameters["functionals"]
            if funs isa AbstractArray
                functionals = Symbol.(funs)
            else
                functionals = [Symbol(funs)]
            end
        end

        if self.parameters["pps"] == "hgh"
            psp_family = "hgh"
            psp_core = :fullcore
        elseif self.parameters["pps"] == "hgh.k"
            psp_family = "hgh"
            psp_core = :semicore
        else
            inputerror_param("pps")
        end

        # Parse smearing and temperature
        temperature = 0.0
        smearing = nothing
        if !isnothing(self.parameters["smearing"])
            psmear = self.parameters["smearing"]
            if lowercase(psmear[1]) == "fermi-dirac"
                smearing = Smearing.FermiDirac()
            elseif lowercase(psmear[1]) == "gaussian"
                smearing = Smearing.Gaussian()
            elseif lowercase(psmear[1]) == "methfessel-paxton"
                smearing = Smearing.MethfesselPaxton(psmear[3])
            else
                inputerror_param("smearing")
            end
            temperature = convert(Float64, psmear[2] / ase_units.Hartree)
        end

        if self.parameters["charge"] != 0.0
            inputerror("Charged systems not supported in DFTK.")
        end

        # TODO This is the place where spin-polarization should be added
        #      once it is in DFTK.

        # Build DFTK atoms
        psploader(symbol) = load_psp(symbol, functional=psp_functional,
                                     family=psp_family, core=psp_core)
        lattice = load_lattice(self.atoms)
        atoms = [ElementPsp(element.symbol, psp=psploader(element.symbol)) => positions
                 for (element, positions) in load_atoms(self.atoms)]

        model_DFT(lattice, atoms, functionals; temperature=temperature,
                  smearing=smearing)
    end

    function get_dftk_basis(self; model=self.get_dftk_model())
        # Build kpoint mesh
        kpts = self.parameters["kpts"]
        kgrid = [1, 1, 1]
        kshift = [0, 0, 0]
        if kpts isa Number
            kgrid = kgrid_size_from_minimal_spacing(model.lattice, kpts * ase_units.Bohr)
        elseif length(kpts) == 3 && all(kpt isa Number for kpt in kpts)
            kgrid = kpts  # Just a plain MP grid
        elseif length(kpts) == 4 && all(kpt isa Number for kpt in kpts[1:3])
            kpts[4] != "gamma" && inputerror("Unknown value to kpts: $kpts")
            kshift = Int.(iseven.([2, 4, 5])) .// 2  # Shift MP grid to always contain Gamma
            kgrid = kpts[1:3]
        elseif kpts isa AbstractArray
            kgrid = nothing
            kshift = nothing
            kcoords = [Vec3(kpt...) for kpt in kpts]
            ksymops = [[(Mat3{Int}(I), Vec3(zeros(3)))] for _ in 1:length(kcoords)]
        end

        # Convert ecut to Hartree
        Ecut = self.parameters["ecut"] / ase_units.Hartree
        if isnothing(kgrid)
            PlaneWaveBasis(model, Ecut, kcoords, ksymops)
        else
            PlaneWaveBasis(model, Ecut, kgrid=kgrid, kshift=kshift)
        end
    end

    function get_dftk_mixing(self; basis=self.get_dftk_basis())
        mixing = basis.model.temperature > 0 ? KerkerMixing(0.7, 1.0) : SimpleMixing(0.7)
        if !isnothing(self.parameters["mixing"])
            if self.parameters["mixing"] isa Tuple
                name = self.parameters["mixing"][1]
                args = ifelse(length(self.parameters["mixing"]) < 2, (),
                              self.parameters["mixing"][2:end])
            else
                name = self.parameters["mixing"]
                args = ()
            end

            valid_types = Dict("KerkerMixing" => KerkerMixing,
                               "SimpleMixing" => SimpleMixing)
            if !haskey(valid_types, name)
                inputerror("A mixing method $name is not known to DFTK.")
            else
                mixing = valid_types[name](args...)
            end
        end
        mixing
    end

    function get_dftk_scfres(self; basis=self.get_dftk_basis())
        mixing = self.get_dftk_mixing(basis=basis)

        extraargs=()
        if !isnothing(self.parameters["nbands"])
            extraargs = (n_bands=self.parameters["nbands"], )
        end

        callback = info -> ()
        if self.parameters["verbose"]
            callback = scf_default_callback
        end

        try
            tol = self.parameters["scftol"]  # TODO Maybe unit conversion here!
            self_consistent_field(basis; tol=tol, callback=callback,
                                  mixing=mixing, extraargs...)
        catch e
            pyraise(calculator.SCFError(string(e)))
        end
    end

    function calculate(self, atoms=nothing, properties=["energy"],
                       system_changes=calculator.all_changes)
        calculator.Calculator.calculate(self, atoms=atoms, properties=properties,
                                        system_changes=system_changes)
        results = self.results  # Load old results

        # Run the SCF and update the energy
        scfres = isnothing(self.scfres) ? self.get_dftk_scfres() : self.scfres
        results["energy"] = sum(values(scfres.energies)) * ase_units.Hartree

        # Compute forces, if requested:
        if "forces" in properties
            # TODO If the calculation fails, ASE expects an
            #      calculator.CalculationFailed exception
            fvalues = hcat((Array(d) for fatom in forces(scfres) for d in fatom)...)'
            results["forces"] = fvalues * (ase_units.Hartree / ase_units.Bohr)
        end

        # Commit results and dump them to disk:
        self.results = results
        self.scfres = scfres
        self.write()
    end
end

end
