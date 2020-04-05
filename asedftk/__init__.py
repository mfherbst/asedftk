from .calculator import DFTK

__all__ = ["DFTK", "install"]

__version__ = "0.1.0"
__license__ = "MIT"
__author__ = ["Michael F. Herbst"]


def register_ase_plugin():
    from ase.calculators.calculator import register_calculator_class

    register_calculator_class("DFTK", DFTK)


def install():
    import julia

    try:
        from julia import Pkg
    except julia.core.UnsupportedPythonError:
        julia.install()

    from julia import Pkg  # noqa: F811

    try:
        from julia import DFTK  # noqa: F401
    except ImportError:
        Pkg.install("DFTK")


register_ase_plugin()
