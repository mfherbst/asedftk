# Setup asedftk from Julia using Julia's conda environment.

# Make a startup config containing the ENV["PYTHON"] = "" in case the user
# has none this far
const CONFIG = abspath(first(DEPOT_PATH), "config", "startup.jl")
if !isfile(CONFIG)
    mkpath(dirname(CONFIG))
    open(CONFIG, "w") do fp
        write(fp, """ENV["PYTHON"] = ""  # Added by setup_asedftk.jl\n""")
    end
elseif !isempty(get(ENV, "PYTHON", ""))
    errorstring = """
    Ensure that your "PYTHON" environment variable is empty to force the usage of
    Julia's internal conda environment for PyCall.jl. One way to enforce this is to add
    the command `ENV["PYTHON"] = ""` to your ~/.julia/config/startup.jl file.
    """
    error(strip(errorstring))
end

println("#\n# Updating Julia dependencies\n#")
ENV["PYTHON"] = ""
ENV["PYTHONPATH"] = ""
import Pkg
Pkg.add(Pkg.PackageSpec(name="Conda", rev="master"))  # For now need Conda from master
Pkg.add(["PyCall", "JSON", "DFTK", "FFTW"])
Pkg.update()
println()

println("#\n# Updating Conda dependencies\n#")
import Conda
Conda.add(["ase"], channel="conda-forge")
Conda.update()
println()

println("#\n# Updating Pip dependencies\n#")
Conda.pip_interop(true)
Conda.pip("install --upgrade", "asedftk")
println()

println("#\nasedftk installation completed in Conda environment $(Conda.ROOTENV)\n#")
println("To use asedftk employ either one of:")
println("    * `conda activate $(Conda.ROOTENV)`\n      followed by `python-jl <script>`")
println("    * directly `$(Conda.ROOTENV)/bin/python-jl <script>`")
