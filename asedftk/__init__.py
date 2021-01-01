import warnings

__version__ = "0.1.7"
__license__ = "MIT"
__author__ = ["Michael F. Herbst"]

# List of all compatible DFTK major.minor versions
COMPATIBLE_DFTK = ["0.1"]


def check_julia():
    """
    Is the 'julia' python package working properly?
    """
    import julia

    try:
        from julia import Main

        julia_compatible = Main.eval('VERSION >= v"1.4.0"')
        if not julia_compatible:
            raise ImportError("Your Julia version is too old. "
                              "At least 1.4.0 required")
        return julia_compatible
    except julia.core.UnsupportedPythonError as e:
        string = ("\n\nIssues between python and Julia. Try to resolve by installing "
                  "required Julia packages using\n"
                  '    python3 -c "import asedftk; asedftk.install()"')
        warnings.warn(str(e) + string)


def dftk_version():
    """
    Get the version of the DFTK package installed Julia-side
    """
    from julia import Main  # noqa: F811, F401
    from julia import Pkg  # noqa: F811, F401

    return Main.eval('''
        string([package.version for (uuid, package) in Pkg.dependencies()
                if package.name == "DFTK"][end])
    ''')


def has_compatible_dftk():
    """
    Do we have a compatible DFTK version installed?
    """
    try:
        from julia import DFTK as jl_dftk  # noqa: F401
    except ImportError:
        return False

    version = dftk_version()
    return any(version.split(".")[:2] == v.split(".") for v in COMPATIBLE_DFTK)


def install(*args, **kwargs):
    import julia

    julia.install(*args, **kwargs)

    try:
        from julia import JSON  # noqa: F401
        from julia import DFTK as jl_dftk  # noqa: F401
    except ImportError:
        from julia import Pkg  # noqa: F811

        Pkg.add("JSON")
        Pkg.add(Pkg.PackageSpec(name="DFTK", version=COMPATIBLE_DFTK[-1]))

from .calculator import DFTK

__all__ = ["install", "dftk_version", "DFTK"]

# if not check_julia():
#     warnings.warn("Julia not found. Try to install Julia requirements "
#                   "using 'asedftk.install()'")
# elif not has_compatible_dftk():
#     warnings.warn("Could not find a compatible DFTK version. Maybe DFTK is not "
#                   "installed on the Julia side or is too old. Before you can "
#                   "use asedftk you have to update it using 'asedftk.install()'")
