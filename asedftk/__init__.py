__version__ = "0.1.0"
__license__ = "MIT"
__author__ = ["Michael F. Herbst"]

# TODO check and ensure a compatible DFTK version is installed.


def has_julia():
    import julia

    try:
        from julia import Main

        return Main.eval("true")
    except julia.core.UnsupportedPythonError:
        return False


__all__ = ["install"]
if has_julia():
    from .calculator import DFTK

    __all__ = ["install", "DFTK"]


def install(*args, **kwargs):
    import julia

    try:
        from julia import Pkg
    except julia.core.UnsupportedPythonError:
        julia.install(*args, **kwargs)

    from julia import Pkg  # noqa: F811

    try:
        from julia import DFTK as jl_dftk  # noqa: F401
    except ImportError:
        Pkg.install("DFTK")
