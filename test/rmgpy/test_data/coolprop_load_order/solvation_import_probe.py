"""
Probe run by ``test/rmgpy/coolpropLoadOrderTest.py`` in a child ``python -m pytest``, so that
``test/conftest.py`` applies to it. Not named ``*Test.py``: the main suite never collects it.
"""


def test_solvation_imports_coolprop():
    import rmgpy.data.solvation as solvation

    assert callable(solvation.PropsSI)
