"""Exceptions shared by every stage."""


class SourceError(RuntimeError):
    """A source file is missing or does not satisfy what the pipeline assumes about it."""


class DataIntegrityError(RuntimeError):
    """Something written or read back does not match what was computed."""


class PipelineStop(RuntimeError):
    """A clean stop (time budget used up, lock held, a dependency missing). Rerun to continue."""
