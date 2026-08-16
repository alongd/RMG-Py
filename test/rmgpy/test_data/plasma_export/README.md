# Golden Cantera exports at the conventional reference temperature

`t0_one_gas.yaml` and `t0_one_surface.yaml` were produced by the Cantera writer
**as it stood at commit `331614fe1`**, i.e. before `rmgpy/yaml_cantera2.py` began
folding the reference temperature `T0` into the exported A-factor.

They exist so that the writer-wide `T0` normalization can be shown to change
nothing for the case RMG actually generates. RMG's convention is `T0 = 1 K`, so
`T0**n == 1` and the division is exactly a division by one; these files are what
turns that from an argument into a check.

Regenerating them is only correct if you first establish, some other way, that
the exported rates have not moved. To regenerate deliberately:

```python
# with rmgpy/yaml_cantera2.py from the reference commit importable as `base`
import plasmaExportTest as t
base.save_cantera_model(t._t0_one_gas_model(t.gas_forms_species.__wrapped__()),
                        'test/rmgpy/test_data/plasma_export/t0_one_gas.yaml')
base.save_cantera_model(t._t0_one_surface_model(t.surface_forms_species.__wrapped__()),
                        'test/rmgpy/test_data/plasma_export/t0_one_surface.yaml',
                        site_density=t._SITE_DENSITY)
```

The comparison normalizes away the `generator:` entry only: it embeds the writer
module's path and the current git commit, neither of which is stable. Everything
else is compared byte for byte.
