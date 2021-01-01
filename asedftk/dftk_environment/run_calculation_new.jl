import JSON
using DFTK
using PyCall
using MPI
ase_units = pyimport("ase.units")

#
# The idea is to use this file to run the caculation and start julia directly from python
# (optionally with mpi processes and threads)
#


function setup_model(params)
    ase_atoms = @pywith pyimport("io").StringIO(params["atoms"]) as f begin
        ase_io.read(f, format="json")
    end
    lattice          = load_lattice(ase_atoms)
    magnetic_moments = load_magnetic_moments(ase_atoms)

    p_psp = params["psp"]
    psploader(symbol) = load_psp(symbol, functional=p_psp["functional"],
                                 family=p_psp["family"], core=Symbol(p_psp["core"]))
    atoms = [ElementPsp(element.symbol, psp=psploader(element.symbol)) => positions
             for (element, positions) in load_atoms(ase_atoms)]

    include_string(Main, params["model"])
end


function run_scf(params, checkpointfile)
    model  = setup_model(params)
    basis  = include_string(Main, params["basis"])
    mixing = include_string(Main, params["mixing"])

    parameters = params["parameters"]
    if !isnothing(parameters["nbands"])
        extraargs = (n_bands=parameters["nbands"], )
    end

    if isfile(checkpointfile)
        oldres    = load_scfres(checkpointfile)
        extraargs = merge(extraargs, (ρspin=oldres.ρspin, ρ=oldres.ρ, ψ=oldres.ψ))
    end

    # TODO Currently taking checkpointing is only implemented if no MPI run
    callback = DFTK.ScfSaveCheckpoints(filename=checkpointfile)
    DFTK.mpi_nprocs() > 1 && (callback = identity)
    if parameters["verbose"]
        callback = callback ∘ DFTK.ScfDefaultCallback()
    end

    self_consistent_field(basis; tol=parameters["scftol"], callback=callback,
                          mixing=mixing, extraargs...)
end

#
# Store results
#
if DFTK.mpi_nprocs() > 1



end


#
# Compute forces
#
# TODO



function main()
    # MPI and threading setup
    DFTK.setup_threading()
    DFTK.mpi_master() || (redirect_stdout(); redirect_stderr())

    @assert endswith(basename, ".json")
    basename, _ = splitext(ARGS[1])

    # Read input data
    if DFTK.mpi_master()
        str = open(fp -> read(fp, String), ARGS[1])
    else
        str = nothing
    end
    MPI.bcast(str, 0, MPI.COMM_WORLD)
    params = open(JSON.parse, ARGS[1], "r")

    scffile = basename * "scfres.jld2"
    if isfile(scffile)
        scfres = load_scfres(scffile)
    else
        checkpointfile = basename * ".checkpoint.jld2"
        scfres = run_scf(params, checkpointfile)

        # TODO Saving scfres in distributed calculations not yet implemented
        if DFTK.mpi_nprocs() > 1
            save_scfres(scffile, scfres)
        end
    end

    # Process results
    results = Dict{String, Any}()
    results["energy"] = scfres.energies.total * ase_units.Hartree

    if "forces" in params["properties"]
        forces = compute_forces_cart(scfres)



            # TODO The next DFTK will have a compute_forces_cart
            #      and the ase_atoms_translation_map function,
            #      which should be used here.

            # DFTK has forces as Hartree over fractional coordinates
            # ASE wants forces as eV / Å
            cart_forces = zeros(eltype(scfres.basis), n_atoms, 3)
            lattice = scfres.basis.model.lattice  # lattice vectors in Bohr
            for (i, atforce) in enumerate(Iterators.flatten(dftk_forces))
                cart_forces[i, :] = lattice \ atforce
            end
            results["forces"] = cart_forces * (ase_units.Hartree / ase_units.Bohr)

    end



    # TODO Write results to json (use the same standard as write and read)


end
