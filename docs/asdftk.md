# Some documentation


## A simple example
Silicon energy and forces:
```python
from asedftk import DFTK
import ase.build

silicon = ase.build.bulk("Si", cubic=True)
silicon.calc = DFTK()

print("Silicon energy: ", silicon.get_potential_energy())
print("Silicon forces: ", silicon.get_forces())
```
