"""Unit conversions for comparing feynlab output with experiment.

Hall conductivities from 2D models come in e^2/h (charge) or e/2pi (= (hbar/e) e^2/h,
spin/orbital).  For a 3D crystal built from layers of spacing d (Angstrom):
    sigma_3D [S/cm] = sigma_2D [e^2/h] * G0 / d,  G0 = e^2/h = 3.874045865e-5 S.
Spin/orbital Hall conductivities are conventionally quoted in (hbar/e) S/cm.
"""
from __future__ import annotations

E2H_SIEMENS = 3.874045865e-5      # e^2/h in S (exact SI, 2019)
HARTREE_EV = 27.211386245988
RYDBERG_EV = 13.605693122994
BOHR_ANGSTROM = 0.529177210903
KB_EV = 8.617333262e-5


def sheet_to_bulk_S_per_cm(sigma_e2h, layer_spacing_angstrom):
    return sigma_e2h * E2H_SIEMENS / (layer_spacing_angstrom * 1e-8)


def bz3d_to_S_per_cm(sigma_e2h_per_angstrom):
    """3D Kubo results integrated as int d^3k/(2pi)^3 with k in 1/A give (e^2/h)/A."""
    return sigma_e2h_per_angstrom * E2H_SIEMENS / 1e-8


def kelvin_to_ev(T):
    return T * KB_EV


def ev_to_kelvin(E):
    return E / KB_EV
