"""Exploratory checks for NASA's Small-Body Database (SBDB) API fields.

The SBDB API (https://ssd-api.jpl.nasa.gov/doc/sbdb.html) can return
physical parameters that NeoWs lacks, e.g., albedo, density, rotation
period, and spectral class. This module documents how we would call the
service and validates the field mapping against a representative sample
response bundled from the public documentation (no live network call is
performed inside the test to respect restricted environments).
"""

from __future__ import annotations

import json
from typing import Any, Dict, Optional

import pytest


SBDB_SAMPLE = {
    "object": {
        "des": "433",
        "fullname": "433 Eros (A898 PA)",
    },
    "orbit": {
        "e": 0.2229687481,
        "a": 1.4581624805,
        "i": 10.8279245056,
        "om": 304.3007558786,
        "w": 178.7492354586,
        "ma": 45.9181555015,
    },
    "phys_par": {
        "H": 10.31,
        "albedo": 0.25,
        "diameter": 16.84,
        "extent": "34.4x11.2x11.2",
        "GM": "0.0004463",
        "mass": "6.69e+15",
        "density": "2670",
        "rot_per": "5.270249",
        "spec": "S",
        "BV": "0.88",
        "UB": "0.39",
    },
}

NEOWS_SAMPLE = {
    "id": "2000433",
    "neo_reference_id": "2000433",
    "name": "433 Eros (A898 PA)",
    "designation": "433",
}


def build_sbdb_request(designation: str, *, include_orbit: bool = True, include_phys: bool = True) -> Dict[str, Any]:
    """Return query parameters for the SBDB API."""

    params: Dict[str, Any] = {"s": designation}
    if include_orbit:
        params["orb"] = "1"
    if include_phys:
        params["phys"] = "1"
    return params


def extract_phys_params(sbdb_payload: Dict[str, Any]) -> Dict[str, Any]:
    """Helper to pull physical parameters from an SBDB response."""

    phys = sbdb_payload.get("phys_par") or {}
    return {
        "H": phys.get("H"),
        "albedo": phys.get("albedo"),
        "diameter_km": phys.get("diameter"),
        "density_kg_m3": _maybe_float(phys.get("density"), scale=1000.0),
        "rotation_period_hr": _maybe_float(phys.get("rot_per")),
        "spectral_class": phys.get("spec"),
    }


def _maybe_float(value: Any, *, scale: Optional[float] = None) -> Optional[float]:
    if value is None:
        return None
    try:
        num = float(value)
    except (TypeError, ValueError):
        return None
    if scale is not None:
        return num * scale
    return num


@pytest.mark.parametrize(
    "designation,expected",
    [
        (
            "433",
            {
                "H": 10.31,
                "albedo": 0.25,
                "diameter_km": 16.84,
                "density_kg_m3": 2.67e6,
                "rotation_period_hr": 5.270249,
                "spectral_class": "S",
            },
        ),
    ],
)
def test_extract_phys_params(designation: str, expected: Dict[str, Any]) -> None:
    # Instead of calling the live API, use the bundled sample payload.
    payload = SBDB_SAMPLE
    assert extract_phys_params(payload) == expected


def test_build_request_flags() -> None:
    assert build_sbdb_request("Apophis") == {"s": "Apophis", "orb": "1", "phys": "1"}
    assert build_sbdb_request("Apophis", include_orbit=False) == {"s": "Apophis", "phys": "1"}
    assert build_sbdb_request("Apophis", include_phys=False) == {"s": "Apophis", "orb": "1"}


def test_cross_reference_designations() -> None:
    sbdb_des = SBDB_SAMPLE["object"]["des"]
    neows_id = NEOWS_SAMPLE["neo_reference_id"]
    neows_designation = NEOWS_SAMPLE.get("designation")
    assert sbdb_des == neows_designation
    # NeoWs reference IDs embed discovery year prefixes; trailing digits match designation.
    assert neows_id.endswith(sbdb_des)


if __name__ == "__main__":
    import requests

    params = build_sbdb_request("99942", include_orbit=True, include_phys=True)
    url = "https://ssd-api.jpl.nasa.gov/sbdb.api"
    response = requests.get(url, params=params, timeout=10)
    response.raise_for_status()
    data = response.json()
    print(json.dumps({
        "designation": data.get("object", {}).get("fullname"),
        "phys_par": extract_phys_params(data),
        "orbit": {k: data.get("orbit", {}).get(k) for k in ("a", "e", "i", "om", "w", "ma")},
    }, indent=2))
