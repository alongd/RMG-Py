"""
A pytest plugin that puts RMG back the way it was before this ticket.

`rmgpy/rmg/model.py` imports `check_quarantine` by name, so replacing that one
binding with a no-op reproduces exactly the pre-ticket code path: the manifest still
loads, the family still knows it is quarantined, and the model admits the reaction
anyway. Run the quarantine tests under it to see which of them are actually load
bearing::

    pytest test/rmgpy/data/kinetics/quarantineTest.py -p no_gate_plugin

Every test that passes under this plugin is one that would have passed before the
gate existed, and so proves nothing on its own.
"""


def pytest_configure(config):
    import rmgpy.rmg.model

    rmgpy.rmg.model.check_quarantine = lambda *args, **kwargs: None
