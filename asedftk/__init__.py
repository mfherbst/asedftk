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

        return Main.eval("true")
    except julia.core.UnsupportedPythonError as e:
        string = ("\n\nIssues between python and Julia. Try to resolve by installing "
                  "required Julia packages using 'asedftk.install()'")
        warnings.warn(e + string)


def dftk_version():
    """
    Get the version of the DFTK package installed Julia-side
    """
    from julia import Main, Pkg  # noqa: F811, F401

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

    return dftk_version() in COMPATIBLE_DFTK


def install(*args, **kwargs):
    import julia

    julia.install(*args, **kwargs)

    from julia import Pkg  # noqa: F811

    try:
        from julia import DFTK as jl_dftk  # noqa: F401
    except ImportError:
        Pkg.install("DFTK@" + COMPATIBLE_DFTK[-1])


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
