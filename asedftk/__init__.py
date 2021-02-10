import io
import os
import json
import shutil
import socket
import datetime
import subprocess
import numpy as np

import ase.io
import ase.build

from ase.calculators.calculator import (CalculationFailed, Calculator,
                                        CalculatorSetupError, Parameters,
                                        all_changes, register_calculator_class)

__version__ = "0.2.5"
__license__ = "MIT"
__author__ = "Michael F. Herbst"
__all__ = ["DFTK", "update", "build_sysimage", "remove_sysimage"]


def environment():
    dftk_environment = os.path.join(os.path.dirname(__file__), "dftk_environment")
    return os.environ.get("ASEDFTK_DFTK_ENVIRONMENT", dftk_environment)


def mpiexecjl():
    # XXX More flexible is to use 'joinpath(DEPOT_PATH[1], "bin", "mpiexecjl")',
    # which is used by MPI to determine the install location.
    return os.path.join(os.path.expanduser("~"), ".julia", "bin", "mpiexecjl")


def sysimagepath():
    return os.path.join(environment(), "sysimage.so")


def julia(*args, n_mpi=1, sysimage=True, **kwargs):
    if sysimage:
        if sysimage is True:
            sysimage = sysimagepath()
        if os.path.isfile(sysimage):
            args = ["--sysimage", sysimage, *args]

    mpiargs = []
    if n_mpi > 1:
        mpiargs = [mpiexecjl(), "--project=" + environment(), "-np", str(n_mpi)]
        args = ["--compiled-modules=no", *args]

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
        return julia("-e", f'VERSION < v"{min_version}" && Sys.exit(1)',
                     sysimage=False)
    except subprocess.CalledProcessError:
        raise RuntimeError(
            f"Julia version below minimal version {min_version}. "
            "Please upgrade julia from https://julialang.org/downloads/."
        )


def build_sysimage():
    """
    Build a custom sysimage for the current julia version and store it for use with
    asedftk. This greatly reduces latency for calculations with asedftk by
    compiling a large part of DFTK into a shared library. This has a few caveats,
    see the asedftk documentation. This feature is currently experimental and may
    not work for you. In case of issues just remove the generated sysimage again
    using `remove_sysimage()`. In this case please file an issue on github to allow
    us to improve support.
    """
    import tempfile

    # Install missing packages, if needed:
    update(always_run=False)

    old = os.getcwd()
    with tempfile.TemporaryDirectory() as tempdir:
        try:
            os.chdir(tempdir)
            precompilefile = os.path.join(os.path.dirname(__file__),
                                          "precompile.jl")
            julia("-e", f"""
                using PackageCompiler;
                create_sysimage([:DFTK, :JLD2, :JSON3],
                sysimage_path="sysimage.so",
                precompile_execution_file="{precompilefile}")
            """, sysimage=False)

            if os.path.isfile(sysimagepath()):
                os.unlink(sysimagepath())
            shutil.move("sysimage.so", sysimagepath())
        finally:
            os.chdir(old)


def remove_sysimage():
    """
    Remove a sysimage created previously with `build_sysimage()`
    """
    if os.path.isfile(sysimagepath()):
        os.unlink(sysimagepath())


def update(always_run=True):
    """
    Update the installed Julia dependencies.
    """
    updated_file = os.path.join(environment(), "last_updated")
    project_file = os.path.join(environment(), "Project.toml")

    # Update file is older than Project file
    needs_update = (
        not os.path.isfile(updated_file)
        or os.stat(updated_file).st_mtime < os.stat(project_file).st_mtime
    )

    if not needs_update and not always_run:
        return

    check_julia_version()
    print("Updating Julia environment ...")
    julia("-e", "import Pkg; Pkg.update(); Pkg.instantiate(); Pkg.precompile()",
          sysimage=False)
    open(updated_file, "w").close()

    if not os.path.isfile(mpiexecjl()):
        julia("-e", "import MPI; MPI.install_mpiexecjl(verbose=true)")
        assert os.path.isfile(mpiexecjl())

    if os.path.isfile(sysimagepath()):
        print("Updating sysimage ...")
        try:
            build_sysimage()
        except subprocess.CalledProcessError as e:
            print(e)
            print("Problems building updated sysimage. You can retry manually "
                  "using 'asedftk.build_sysimage()'.")


def run_calculation(properties, inputfile, n_threads=1, n_mpi=1):
    check_julia_version()
    update(always_run=False)
    script = os.path.join(os.path.dirname(__file__), "run_calculation.jl")
    logfile = os.path.splitext(inputfile)[0] + ".log"

    try:
        with open(logfile, "a") as fp:
            fp.write("#\n")
            fp.write(f"#--  {socket.gethostname()}  {datetime.datetime.now()}\n")
            fp.write("#\n")
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
            "n_threads": 1,
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
        self.scfres = saved_dict.get("scfres", None)
        self.results = saved_dict["results"]

        # Some results need to be numpy arrays:
        for key in ("forces", ):
            self.results[key] = np.array(self.results[key])

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
        inputfile = os.path.abspath(self.label + ".json")
        self.write(label=self.label)

        # Run DFTK
        n_threads = self.parameters["n_threads"]
        n_mpi = self.parameters["n_mpi"]
        run_calculation(properties, inputfile, n_threads=n_threads,
                        n_mpi=n_mpi)

        # Read results
        self.read(self.label)


# ... and register it with ASE
register_calculator_class("DFTK", DFTK)
