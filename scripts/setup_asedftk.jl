# Setup asedftk from Julia using Julia's conda environment.

# Make a startup config containing the ENV["PYTHON"] = "" in case the user
# has none this far
const CONFIG = abspath(first(DEPOT_PATH), "config")
if !isfile(join(CONFIG, "startup.jl"))
    mkpath(CONFIG)
    open("startup.jl", "w") do fp
        write(fp, raw"ENV["PYTHON"] = ""  # Added by setup_asedftk.jl\n")
    end
elseif !isempty(get(ENV, "PYTHON", ""))
    errorstring = """
    Ensure that your "PYTHON" environment variable is empty to force the usage of
    Julia's internal conda environment for PyCall.jl. One way to enforce this is to add
    the command `ENV["PYTHON"] = ""` to your ~/.julia/config/startup.jl file.
    """
    error(strip(errorstring))
end

# Setup Conda and Julia packages
ENV["PYTHON"] = ""
ENV["PYTHONPATH"] = ""
import Pkg
Pkg.add(Pkg.PackageSpec(name="Conda", rev="master"))  # For now need Conda from master
Pkg.add(["PyCall", "JSON", "DFTK", "FFTW"])
Pkg.update()

# Add/update conda dependencies
import Conda
Conda.add(["ase"], channel="conda-forge")
Conda.update()

# Add/update pip dependencies
if PyCall.ispynull(pyimport_e("asedftk"))
    Conda.pip_interop(true)
    Conda.pip("install", "asedftk")
else
    Conda.pip("install --upgrade --no-deps", "asedftk")
end

println("asedftk installation completed in Conda environment $(Conda.ROOTENV)")
println("To use asedftk employ either one of:")
println("    * `conda activate $(Conda.ROOTENV)`; followed by `python-jl <script>`")
println("    * directly `$(Conda.ROOTENV)/bin/python-jl <script>`")
