"""
Public material-response helpers re-exported for the API surface.

Workflow role
-------------
The package keeps the low-level constitutive routines split across invariant and
plasticity modules, but the public API exposes the functions from a single
location. This page is therefore the best entry point for readers who want to
see which stress-update utilities are intended for direct use from outside the
element kernels.
"""

from .invariants import devstress, eqstress
from .plasticity import dyieldvm, stressdp, stressvm, yieldvm

devstres = devstress

__all__ = [
    "devstress",
    "devstres",
    "dyieldvm",
    "eqstress",
    "stressdp",
    "stressvm",
    "yieldvm",
]
