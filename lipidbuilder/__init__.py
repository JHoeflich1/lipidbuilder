__all__ = ["LipidSystemBuilder"]


def __getattr__(name):
    if name == "LipidSystemBuilder":
        from .builder import LipidSystemBuilder

        return LipidSystemBuilder
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
