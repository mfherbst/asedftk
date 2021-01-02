import os

from asedftk import julia, DFTK
import ase
import ase.build


def dump_data():
    # Dump atoms data
    cell = [[3.21, 0.0, 0.0], [-1.605, 2.7799415461480477, 0.0], [0.0, 0.0, 5.21304]]
    positions = [[0.0, 0.0, 0.0], [0.0, 1.85329, 2.60652]]
    atoms = ase.Atoms(symbols='Mg2', pbc=True, cell=cell, positions=positions)
    calc = DFTK()
    calc.atoms = atoms
    calc.write(label="mg")

    calc = DFTK()
    calc.atoms = ase.build.bulk("Si")
    calc.write(label="si")


def test_juliatests():
    dump_data()

    # Run tests
    thisdir = os.path.dirname(__file__)
    julia(os.path.join(thisdir, "runtests.jl"))
