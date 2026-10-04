"""Variant identity. No pathogenicity convenience method."""

from __future__ import annotations

from typing import Any, Mapping

from modules.variant_core import build_variant, parse_variant_input
from modules.variant_explorer import explore_variant as explore_variant_layers

__all__ = ("build_variant", "explore_variant", "identify_variant", "parse_variant_input")


def identify_variant(text: str, *, assembly: str) -> dict[str, Any]:
    """Normalize variant identity. Assembly is mandatory. Effect is not computed."""
    return build_variant(text=text, assembly=assembly)


def explore_variant(
    *,
    variant: Mapping[str, Any],
    email: str = "",
    enable_vep: bool = True,
    enable_clinvar: bool = True,
    enable_domains: bool = True,
    vep_urlopen_fn=None,
    clinvar_urlopen_fn=None,
    domain_urlopen_fn=None,
) -> dict[str, Any]:
    """Remote annotation layers. Does not diagnose pathogenicity."""
    return explore_variant_layers(
        variant=variant,
        email=email,
        enable_vep=enable_vep,
        enable_clinvar=enable_clinvar,
        enable_domains=enable_domains,
        vep_urlopen_fn=vep_urlopen_fn,
        clinvar_urlopen_fn=clinvar_urlopen_fn,
        domain_urlopen_fn=domain_urlopen_fn,
    )
