import ase
import asedftk


def test_registered():
    assert ase.calculators.calculator.get_calculator_class("DFTK") == asedftk.DFTK
