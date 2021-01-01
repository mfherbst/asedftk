import ase.io
import io
import os
import json

from ase.calculators.calculator import Calculator, CalculatorSetupError, Parameters
from ase.calculators.calculator import register_calculator_class, all_changes

class DFTK(Calculator):
    implemented_properties = ["energy", "forces"]

    def __init__(self, atoms=None, label="dftk", **kwargs):
        self.default_parameters = {
            "xc": "LDA",
            "kpts": 3.5,
            "smearing": None,
            "nbands": None,
            "charge": 0.0,
            "functionals": None,
            "pps": "hgh",
            "scftol": 1e-5,
            "ecut": 400,
            "verbose": False,
            "mixing": None,
        }
        super().__init__(label=label, **kwargs)
        self.scfres = None

    def reset(self):
        # Reset internal state by purging all intermediates
        super().reset()
        if self.scfres and os.path.isfile(self.scfres):
            os.unlink(self.scfres)
            self.scfres = None

    def read(self, label):
        super().read(label)

        with open(label + ".json", "r") as fp:
            saved_dict = json.load(fp)
        self.parameters = Parameters(saved_dict["parameters"])
        self.results = saved_dict["results"]
        self.scfres = saved_dict.get("scfres", None)

        with io.StringIO(saved_dict["atoms"]) as fp:
            self.atoms = ase.io.read(fp, format="json")

    def write(self, label=None):
        if label is None:
            label = self.label
        if self.atoms is None:
            raise CalculatorSetupError("An Atoms object must be present "
                                       "in the calculator to use write")

        with io.StringIO() as fp:
            ase.io.write(fp, self.atoms.copy(), format="json")
            atoms_json = fp.getvalue()

        with open(label + ".json", "w") as fp:
            save_dict = {
                "parameters": self.parameters,
                "results": self.results,
                "atoms": atoms_json,
                "scfres": self.scfres,
            }
            json.dump(save_dict, fp)

    def set(self, **kwargs):
        changed_parameters = super().set(**kwargs)
        if changed_parameters:
            # Currently reset on all parameter changes
            # TODO Improve
            self.reset()
        return changed_parameters

    def check_state(self, atoms):
        system_changes = super().check_state(atoms)
        if "pbc" in system_changes:
            # Ignore boundary conditions (we always use PBCs)
            system_changes.remove("pbc")
        return system_changes

    def calculate(self, atoms=None, properties=["energy"], system_changes=all_changes):
        super().calculate(atoms=atoms, properties=properties, system_changes=system_changes)

        # Write input file
        inputfile = self.label + ".json"
        self.write()

        # Run DFTK
        commandline = "julia run_dftk.jl " + " ".join(properties) + " " + inputfile
        # TODO Run commandline
        # TODO mpirun, -nt --project

        # Read results
        self.read(self.label)


# ... and register it with ASE
register_calculator_class("DFTK", DFTK)
