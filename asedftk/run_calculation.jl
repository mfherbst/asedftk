using JSON
using DFTK
using JLD2
using MPI
using Unitful
using UnitfulAtomic

function setup()
    if DFTK.mpi_nprocs() > 1
        DFTK.disable_threading()
    else
        DFTK.setup_threading()
    end
    DFTK.mpi_master() || (redirect_stdout(); redirect_stderr())
end

inputerror(parameters, p) = error("Unknown value to parameter $p: $(get(parameters, p, ""))")


function get_dftk_model(parameters, extra)
    # Parse psp and DFT functional parameters
    functionals = [:lda_xc_teter93]
    if lowercase(parameters["xc"]) == "lda"
        psp_functional = "lda"
        functionals = [:lda_xc_teter93]
    elseif lowercase(parameters["xc"]) == "pbe"
        psp_functional = "pbe"
        functionals = [:gga_x_pbe, :gga_c_pbe]
    else
        inputerror(parameters, "xc")
    end

    if !isnothing(parameters["functionals"])
        funs = parameters["functionals"]
        if funs isa AbstractArray
            functionals = Symbol.(funs)
        else
            functionals = [Symbol(funs)]
        end
    end

    if parameters["pps"] == "hgh"
        psp_family = "hgh"
        psp_core = :fullcore
    elseif parameters["pps"] == "hgh.k"
        psp_family = "hgh"
        psp_core = :semicore
    else
        inputerror(parameters, "pps")
    end

    # Parse smearing and temperature
    temperature = 0.0u"eV"
    smearing = nothing
    if !isnothing(parameters["smearing"])
        psmear = parameters["smearing"]
        if lowercase(psmear[1]) == "fermi-dirac"
            smearing = Smearing.FermiDirac()
        elseif lowercase(psmear[1]) == "gaussian"
            smearing = Smearing.Gaussian()
        elseif lowercase(psmear[1]) == "methfessel-paxton"
            smearing = Smearing.MethfesselPaxton(psmear[3])
        else
            inputerror(parameters, "smearing")
        end
        temperature = convert(Float64, psmear[2]) * u"eV"
    end

    if parameters["charge"] != 0.0
        error("Charged systems not supported in DFTK.")
    end

    # TODO This is the place where spin-polarization should be added
    #      once it is in DFTK.

    # Build DFTK atoms
    psploader(symbol) = load_psp(symbol, functional=psp_functional,
                                 family=psp_family, core=psp_core)
    atoms = [ElementPsp(element.symbol, psp=psploader(element.symbol)) => positions
             for (element, positions) in extra["atoms"]]

    model_DFT(extra["lattice"], atoms, functionals; temperature=temperature,
              smearing=smearing)
end


function get_dftk_basis(parameters, extra; model=get_dftk_model(parameters, extra))
    # Build kpoint mesh
    kpts = parameters["kpts"]
    kgrid = [1, 1, 1]
    kshift = [0, 0, 0]
    if kpts isa Number
        # DFTK uses k-point spacing whereas ASE uses k-point density
        kgrid = kgrid_size_from_minimal_spacing(model.lattice, 1 / kpts / u"Å")
    elseif length(kpts) == 3 && all(kpt isa Number for kpt in kpts)
        kgrid = kpts  # Just a plain MP grid
    elseif length(kpts) == 4 && all(kpt isa Number for kpt in kpts[1:3])
        kpts[4] != "gamma" && error("Unknown value to kpts: $kpts")
        kgrid = kpts[1:3]
        kshift = Int.(iseven.(kgrid)) .// 2  # Shift MP grid to always contain Gamma
    elseif kpts isa AbstractArray
        kgrid = nothing
        kshift = nothing
        kcoords = [Vec3(kpt...) for kpt in kpts]
        ksymops = [[(Mat3{Int}(I), Vec3(zeros(3)))] for _ in 1:length(kcoords)]
    end

    # Convert ecut to Hartree
    Ecut = parameters["ecut"] * u"eV"
    if isnothing(kgrid)
        PlaneWaveBasis(model, austrip(Ecut), kcoords, ksymops)
    else
        PlaneWaveBasis(model, Ecut, kgrid=kgrid, kshift=kshift)
    end
end


function get_dftk_mixing(parameters, extra; basis=get_dftk_basis(parameters, extra))
    if isnothing(parameters["mixing"])
        basis.model.temperature > 0 ? KerkerMixing() : SimpleMixing(α=0.8)
    else
        include_string(Main, parameters["mixing"], @__FILE__)
    end
end


function get_dftk_scfres(parameters, extra; basis=get_dftk_basis(parameters, extra))
    mixing = get_dftk_mixing(parameters, extra; basis=basis)

    extraargs = ()
    if !isnothing(parameters["nbands"])
        extraargs = (n_bands=parameters["nbands"], )
    end

    callback = DFTK.ScfSaveCheckpoints(extra["checkpointfile"]) ∘ DFTK.ScfDefaultCallback()
    self_consistent_field(basis; tol=parameters["scftol"], callback=callback,
                          mixing=mixing, extraargs...)
end

function parse_json_array(data::Dict)
    # Simplistic parser for ndarrays as they are stored in ASE json files
    shape, type, values = data["__ndarray__"]
    shape = reverse(Int.(shape))  # row-major -> column-major
    if type == "float64"
        reshape(Float64.(values), shape...)
    elseif type == "int64"
        reshape(Int64.(values), shape...)
    else
        error("Type not implemented: $type")
    end
end


function load_state(file)
    endswith(file, ".json") || error("State file should end in .json")

    # Read input data
    if DFTK.mpi_master()
        str = open(fp -> read(fp, String), file)
    else
        str = nothing
    end
    MPI.bcast(str, 0, MPI.COMM_WORLD)
    res = JSON.parse(str)

    # Parse atoms json
    atoms_json = JSON.parse(res["atoms"])
    if length(atoms_json["ids"]) != 1
        @warn "Only parsing last atoms object in json. Ignoring all others"
    end
    atoms_id = string(atoms_json["ids"][end])

    lattice_Ang   = parse_json_array(atoms_json[atoms_id]["cell"]["array"])
    positions_Ang = parse_json_array(atoms_json[atoms_id]["positions"])
    numbers       = parse_json_array(atoms_json[atoms_id]["numbers"])

    atoms = []
    positions_fractional = Array.(eachcol(lattice_Ang \ positions_Ang))
    for (i, pos) in enumerate(positions_fractional)
        id = findfirst(el_pos -> el_pos[1].Z == numbers[i], atoms)
        if isnothing(id)
            push!(atoms, ElementCoulomb(numbers[i]) => [pos])
        else
            push!(atoms[id][2], pos)
        end
    end

    res["extra"] = Dict{String, Any}()
    res["extra"]["lattice"]    = lattice_Ang * austrip(1u"Å")
    res["extra"]["atoms"]      = atoms
    res
end


function save_state(file, state)
    endswith(file, ".json") || error("State file should end in .json")

    # Transpose arrays inside results, such that they are written
    # in a predictable order to disk
    save_results = Dict{String, Any}()
    for (key, val) in pairs(state["results"])
        if val isa AbstractArray
            if ndims(val) == 2
                save_results[key] = Vector.(eachrow(val))
                continue
            elseif ndims(val) > 2
                error("Arrays of dimension larger 2 not supported")
            end
        end
        save_results[key] = val
    end

    save_dict = Dict(
        "parameters" => state["parameters"],
        "results"    => save_results,
        "atoms"      => state["atoms"],
        "scfres"     => state["scfres"],
    )
    if DFTK.mpi_master()
        open(file, "w") do fp
            JSON.print(fp, save_dict)
        end
    end
    save_dict
end


function run_calculation(properties::AbstractArray, statefile::AbstractString)
    state = load_state(statefile)
    @assert endswith(state["scfres"], ".scfres.jld2")
    prefix = state["scfres"][1:end-12]
    state["extra"]["checkpointfile"] = prefix * ".checkpoint.jld2"

    if !isfile(state["scfres"])
        save_scfres(state["scfres"], get_dftk_scfres(state["parameters"], state["extra"]))
    end
    scfres = load_scfres(state["scfres"])

    state["results"]["energy"] = ustrip(auconvert(u"eV", scfres.energies.total))
    if "forces" in properties
        # TODO If the calculation fails, ASE expects an
        #      calculator.CalculationFailed exception
        forces = compute_forces_cart(scfres)

        # DFTK has forces as Hartree over fractional coordinates
        # ASE wants forces as eV / Å
        n_atoms = sum(length, forces)
        cart_forces = zeros(eltype(scfres.basis), n_atoms, 3)
        for (i, atforce) in enumerate(Iterators.flatten(forces))
            cart_forces[i, :] = @. ustrip(auconvert(u"eV/Å", atforce))
        end
        state["results"]["forces"] = cart_forces
    end
    save_state(statefile, state)
end


function main()
    setup()
    run_calculation(ARGS[1:end-1], ARGS[end])
end

(abspath(PROGRAM_FILE) == @__FILE__) && main()
