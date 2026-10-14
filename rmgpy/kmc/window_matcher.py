"""Window matcher (red stub): templates + MeltState -> transitions."""

RESONANCE_FORM_LIMIT = 64


class WindowMatcherError(ValueError):
    pass


class TemplateError(WindowMatcherError):
    pass


class ResonanceLimitError(WindowMatcherError):
    pass


class ConjugatedRingError(ResonanceLimitError):
    pass


class StaleTransitionError(WindowMatcherError):
    pass


def _todo(*args, **kwargs):
    raise NotImplementedError("window matcher not implemented yet")


WindowAtom = BondEdit = AtomEdit = PathPredicate = Template = PairSampler = _todo
apply_transition = component_forms = form_rows = match = product_rows = _todo
