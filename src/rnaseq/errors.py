"""Domain-specific exceptions for expected project failures."""


class ProjectConfigError(ValueError):
    """Raised when project configuration cannot be parsed or validated."""


class ProjectCreationError(ValueError):
    """Raised when a project cannot be safely created."""


class ExecutionPreflightError(RuntimeError):
    """Raised when a project cannot safely begin upstream execution."""


class UpstreamExecutionError(RuntimeError):
    """Raised when the Nextflow/nf-core upstream run does not complete."""


class DownstreamExecutionError(RuntimeError):
    """Raised when the bounded L1 expression-QC stage cannot complete."""
