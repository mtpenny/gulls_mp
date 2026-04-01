"""Helpers for running gulls smoke tests."""


def main(*args, **kwargs):
    """Lazy wrapper to avoid importing plotting-heavy modules on package import."""
    from .runner import main as _main

    return _main(*args, **kwargs)


__all__ = ["main"]
