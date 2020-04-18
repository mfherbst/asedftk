# ASE-compatible calculator using DFTK
[![Build Status](https://api.travis-ci.com/mfherbst/asedftk.svg?branch=master)](https://travis-ci.com/mfherbst/asedftk)

At least julia 1.3, 1.4 recommended

Setup:
```
pip install asedftk

import asedftk
asedftk.install()
```

Usage:
```
from asedftk import DFTK

...

atoms.calc = DFTK()

atoms.get_potential_energy()
```

Tips and tricks:
- Threading:
  ```
  JULIA_NUM_THREADS=$NCPUS
  ```
- `python-jl` instead of actual python


See [documentation](docs/asedftk.md).
