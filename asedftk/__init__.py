__version__ = "0.2.2"
__license__ = "MIT"
__author__ = ["Michael F. Herbst"]

import io
import os
import json
import datetime
import subprocess

import ase.io

from ase.calculators.calculator import (CalculationFailed, Calculator,
                                        CalculatorSetupError, Parameters,
                                        all_changes, register_calculator_class)

__all__ = ["DFTK", "update"]


def environment():
    dftk_environment = os.path.join(os.path.dirname(__file__), "dftk_environment")
    return os.environ.get("ASEDFTK_DFTK_ENVIRONMENT", dftk_environment)


def julia(*args, n_mpi=1, **kwargs):
    mpiargs = []
    if n_mpi > 1:
        raise NotImplementedError("MPI not yet implemented. Please select n_mpi=1.")
        mpiargs = [get_mpiexecjl(), "--project=" + environment(), "-np", str(n_mpi)]

    julia_exe = os.environ.get("JULIA", "julia")
    try:
        args = mpiargs + [julia_exe, "--project=" + environment(),
                          "--startup-file=no", *args]
        return subprocess.check_call(args, **kwargs)
    except FileNotFoundError:
        raise RuntimeError(
            "Julia not found. Please check that Julia is installed and "
            "added to your PATH variable or otherwise point the JULIA "
            "environment variable to the full path of your julia binary."
        )


def check_julia_version(min_version="1.5.0"):
    try:
        return julia("-e", f'VERSION < v"{min_version}" && Sys.exit(1)')
    except subprocess.CalledProcessError:
        raise RuntimeError(
            f"Julia version below minimal version {min_version}. "
            "Please upgrade."
        )


def get_mpiexecjl():
    # XXX More flexible is to use 'joinpath(DEPOT_PATH[1], "bin", "mpiexecjl")',
    # which is used by MPI to determine the install location.
    return os.path.join(os.path.expanduser("~"), ".julia", "bin", "mpiexecjl")


def update(always_run=True):
    """
    Update the installed Julia dependencies.
    """
    updated_file = os.path.join(environment(), "last_updated")
    project_file = os.path.join(environment(), "Project.toml")

    # Update file is older than Project file
    needs_update = (not os.path.isfile(updated_file)
                    or os.stat(updated_file).st_mtime < os.stat(project_file).st_mtime)

    if not needs_update and not always_run:
        return

    check_julia_version()
    julia("-e", "import Pkg; Pkg.instantiate(); Pkg.update(); Pkg.precompile()")
    open(updated_file, "w").close()

    if not os.path.isfile(get_mpiexecjl()):
        julia("-e", "import MPI; MPI.install_mpiexecjl(verbose=true)")
        assert os.path.isfile(get_mpiexecjl())


def run_calculation(properties, inputfile, n_threads=1, n_mpi=1):
    check_julia_version()
    update(always_run=False)
    script = os.path.join(os.path.dirname(__file__), "run_calculation.jl")
    logfile = os.path.splitext(inputfile)[0] + ".log"
    if n_threads is None:
        try:
            n_threads = int(os.environ.get("JULIA_NUM_THREADS", "1"))
        except ValueError:
            n_threads = 1

    try:
        with open(logfile, "a") as fp:
            fp.write(f"#\n#--  {datetime.datetime.now()}\n#\n")
            fp.flush()
            julia(*["-t", str(n_threads), script, *properties, inputfile],
                  n_mpi=n_mpi, stderr=subprocess.STDOUT, stdout=fp)
            fp.flush()
            fp.write("\n")
    except subprocess.CalledProcessError:
        msg = f"DFTK calculation failed. See logfile {logfile} for details."
        if "CI" in os.environ:
            with open(logfile, "r") as fp:
                msg += "\n\n" + fp.read()
        raise CalculationFailed(msg)


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
            "mixing": None,
            "n_mpi": 1,
            "n_threads": None,
        }
        self.scfres = label + ".scfres.jld2"
        super().__init__(label=label, atoms=atoms, **kwargs)

    def reset(self):
        # Reset internal state by purging all intermediates
        super().reset()
        if os.path.isfile(self.scfres):
            os.unlink(self.scfres)

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

    def calculate(self, atoms=None, properties=["energy"],
                  system_changes=all_changes):
        super().calculate(atoms=atoms, properties=properties,
                          system_changes=system_changes)

        # Write input file
        inputfile = self.label + ".json"
        self.write()

        # Run DFTK
        n_threads = self.parameters["n_threads"]
        n_mpi = self.parameters["n_mpi"]
        run_calculation(properties, inputfile, n_threads=n_threads,
                        n_mpi=n_mpi)

        # Read results
        self.read(self.label)


# ... and register it with ASE
register_calculator_class("DFTK", DFTK)
