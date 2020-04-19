import warnings

__version__ = "0.1.0"
__license__ = "MIT"
__author__ = ["Michael F. Herbst"]

# List of all compatible DFTK versions
COMPATIBLE_DFTK = ["0.0.6"]


def check_julia():
    """
    Is the 'julia' python package working properly?
    """
    import julia

    try:
        from julia import Main

        julia_compatible = Main.eval('VERSION >= v"1.3.0"')
        if not julia_compatible:
            raise ImportError("Your Julia version is too old. "
                              "At least 1.3.0 required")
        return julia_compatible
    except julia.core.UnsupportedPythonError as e:
        string = ("\n\nIssues between python and Julia. Try to resolve by installing "
                  "required Julia packages using\n"
                  '    python -m "import asedftk; asedftk.install()"')
        warnings.warn(str(e) + string)


def dftk_version():
    """
    Get the version of the DFTK package installed Julia-side
    """
    from julia import Main, Pkg  # noqa: F811, F401

    if Main.eval('VERSION >= v"1.4.0"'):
        return Main.eval('''
            string([package.version for (uuid, package) in Pkg.dependencies()
                    if package.name == "DFTK"][end])
        ''')
    else:
        warnings.warn("Cannot determine DFTK version if using Julia < 1.4.0. "
                      "Hoping for the best ...")
        return COMPATIBLE_DFTK[-1]


def has_compatible_dftk():
    """
    Do we have a compatible DFTK version installed?
    """
    try:
        from julia import DFTK as jl_dftk  # noqa: F401
    except ImportError:
        return False

    return dftk_version() in COMPATIBLE_DFTK


def install(*args, **kwargs):
    import julia

    julia.install(*args, **kwargs)

    try:
        from julia import DFTK as jl_dftk, JSON  # noqa: F401
    except ImportError:
        from julia import Pkg  # noqa: F811

        Pkg.add("JSON")
        Pkg.Registry.add(Pkg.RegistrySpec(
            url="https://github.com/JuliaMolSim/MolSim.git"))
        Pkg.add(Pkg.PackageSpec(name="DFTK", version=COMPATIBLE_DFTK[-1]))


__all__ = ["install", "dftk_version"]
if not check_julia():
    warnings.warn("Julia not found. Try to install Julia requirements "
                  "using 'asedftk.install()'")
elif not has_compatible_dftk():
    warnings.warn("Could not find a compatible DFTK version. Maybe DFTK is not "
                  "installed on the Julia side or is too old. Before you can "
                  "use asedftk you have to update it using 'asedftk.install()'")
else:
    from .calculator import DFTK

    __all__ = ["install", "dftk_version", "DFTK"]
