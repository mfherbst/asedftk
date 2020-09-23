import os

from julia import Main
from ase.calculators.calculator import register_calculator_class

# Get calculator from Julia side ...
script_dir = os.path.dirname(os.path.realpath(__file__))
DFTK = Main.include(os.path.join(script_dir, "calculator.jl"))

# ... and register it with ASE
register_calculator_class("DFTK", DFTK)
