# ASE-compatible calculator using DFTK

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
