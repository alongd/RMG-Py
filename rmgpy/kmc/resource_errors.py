"""Resource errors are fatal (I-265).

A MemoryError or a RecursionError says that the process ran out of memory or
stack, not what the chemistry is.  A handler that turns failures into
refusals, gaps, error certificates or error rows made such an error a verdict
that depends on the machine: two emission shards of one artifact certified a
record differently, republished different ids, and sliced different record
lists.  Every broad handler on the compile, certification, emission and
assembly paths re-raises them first (``reraise_resource``), so the process
fails with a non-zero exit instead.
"""

RESOURCE_ERRORS = (MemoryError, RecursionError)


def reraise_resource(error: BaseException) -> None:
    """Re-raise ``error`` when it is a resource error; return otherwise."""
    if isinstance(error, RESOURCE_ERRORS):
        raise error
