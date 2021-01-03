# Installation of asedftk

TODO Update

## Option 1: Using Julia's conda environment (Recommended)

This first option is the recommended way to get started with asedftk.
It is best suited for users that do not yet have experience with Julia
and do not yet have a working Julia installation.

1. Install Julia e.g. by [downloading the binary](https://julialang.org/downloads).
   The use of at least Julia **1.5** is required.
1. Install asedftk using the [setup_asedftk.jl](https://raw.githubusercontent.com/mfherbst/asedftk/master/scripts/setup_asedftk.jl) script.
   For this [download the setup_asedftk.jl script](https://raw.githubusercontent.com/mfherbst/asedftk/master/scripts/setup_asedftk.jl) and run it with Julia,
   i.e. execute `/path/to/bin/julia setup_asedftk.jl`.
   This automatically installs required Julia dependencies
   such as [DFTK.jl](https://dftk.org) and Python dependencies
   such as [PyJulia](https://pypi.org/project/julia/),
   which is a package to allow
   Julia and Python codes to interoperate with each other.
   The python dependencies are installed inside the conda environment
   managed by [Conda.jl](https://github.com/JuliaPy/Conda.jl/).
   If you have already used your Julia installation beforehand,
   this step might have side effects, see the note below.
1. To use asedftk first activate Julia's conda environment using
   `conda activate ~/.julia/conda/3`, then use the `python-jl` wrapper
   whenever you want to run a python script making use of asedftk.
   I.e. if you have written a calculation script `script.py` you
   need to start it as `python-jl script.py` in order to be able to properly
   use asedftk.

To update asedftk at a later point simply run same script again,
call `/path/to/julia setup_asedftk.jl`.

**Note:** This installation has side effects if you are already using Conda.jl
or [PyCall.jl](https://github.com/JuliaPy/PyCall.jl)
from your Julia installation. More precisely it will rebuild
PyCall.jl in order to use the conda environment shipped by Conda.jl
and will install some packages in the latter environment.
This implies that PyCall.jl and the `pyimport` functions will not use the
system or user package repositories (the ones managed with `pip`), but
use the conda environment from Conda.jl. If this is not what you want
see option 2.


## Option 2: Using an external python environment

This installation assumes that you already have Julia (at least **1.4**) installed
and the [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) package
uses an external python environment, where the python packages are
managed using `pip`.

1. Install asedftk from [PyPi](https://pypi.org/project/asedftk) into your python site packages:
   ```
   pip install asedftk
   ```
   This automatically installs the [PyJulia](https://pypi.org/project/julia/) package,
   which allows Julia and Python codes to interoperate with each other.
1. Install the Julia dependencies of asedftk:
   ```
   python3 -c "import asedftk; asedftk.install()"
   ```
1. That's it, you're all set. But **please note**:
   Due to [some limitations](https://pyjulia.readthedocs.io/en/stable/troubleshooting.html#your-python-interpreter-is-statically-linked-to-libpython)
   in some Linux distros like Debian or Ubuntu
   you might need to run your Python scripts
   with the `python-jl` wrapper if you want to use asedftk in them.
   I.e. if you have written a calculation script `script.py` you
   might need to start it as `python-jl script.py`
   in order to be able to use asedftk.
