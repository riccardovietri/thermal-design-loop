"""Coolant properties at a fixed temperature, from CoolProp (Melinder incompressible data for glycols)."""
from CoolProp.CoolProp import PropsSI

from .types import Fluid

# PG25 is 25% propylene glycol by volume, about 26% by mass.
_COOLPROP_NAMES = {"water": "Water", "PG25": "INCOMP::MPG[0.26]"}


def fluid_props(name: str, temp_c: float) -> Fluid:
    cp_name = _COOLPROP_NAMES[name]
    t, p = temp_c + 273.15, 2.0e5
    rho, cp, k, mu = (PropsSI(x, "T", t, "P", p, cp_name) for x in ("D", "C", "L", "V"))
    return Fluid(name, temp_c, rho, cp, k, mu)
