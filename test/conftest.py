import multiprocessing

multiprocessing.set_start_method('fork')

# rmgpy must be imported before pybel: rmgpy/__init__.py loads CoolProp before any openbabel
# format plugin, which pybel loads at import (see the comment there).
import rmgpy  # noqa: F401

from openbabel import pybel

pybel.ob.obErrorLog.SetOutputLevel(0)
