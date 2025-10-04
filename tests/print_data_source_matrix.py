"""Utility script to inspect white-paper parameters against NeoWs/SBDB availability.

Run with:

    python3 tests/print_data_source_matrix.py

This prints each PAIR-model parameter, whether it appears in the NeoWs
sample payload, the SBDB sample payload, or neither. Live API calls are
optional; when `LIVE=1` is exported and network is permitted, the script
will attempt to retrieve the same fields for a specific asteroid.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Dict, Optional

import requests


# NeoWs sample payload (subset) paralleling what fetch_today_neos would return.
NEOWS_SAMPLE = {
    "id": "2000433",
    "neo_reference_id": "2000433",
    "name": "433 Eros (A898 PA)",
    "designation": "433",
    "absolute_magnitude_h": 10.4,
    "estimated_diameter": {
        "meters": {
            "estimated_diameter_min": 15580.0,
            "estimated_diameter_max": 17460.0,
        }
    },
    "is_potentially_hazardous_asteroid": False,
    "close_approach_data": [
        {
            "relative_velocity": {
                "kilometers_per_second": "5.576",
            },
            "miss_distance": {
                "kilometers": "643234.0",
            },
            "orbiting_body": "Earth",
        }
    ],
}


SBDB_SAMPLE = {
    "object": {
        "des": "433",
        "fullname": "433 Eros (A898 PA)",
    },
    "orbit": {
        "a": 1.4581624805,
        "e": 0.2229687481,
        "i": 10.8279245056,
        "om": 304.3007558786,
        "w": 178.7492354586,
        "ma": 45.9181555015,
    },
    "phys_par": {
        "H": 10.31,
        "albedo": 0.25,
        "diameter": 16.84,
        "density": "2670",
        "rot_per": "5.270249",
        "spec": "S",
    },
}


@dataclass
class ParameterStatus:
    name: str
    source: str
    available_neows: bool
    available_sbdb: bool
    notes: str


PARAMETERS = [
    ParameterStatus(
        name="Absolute magnitude (H)",
        source="PAIR §2.1",
        available_neows="absolute_magnitude_h" in NEOWS_SAMPLE,
        available_sbdb="H" in SBDB_SAMPLE.get("phys_par", {}),
        notes="NeoWs + SBDB",
    ),
    ParameterStatus(
        name="Diameter",
        source="PAIR §2.1",
        available_neows="estimated_diameter" in NEOWS_SAMPLE,
        available_sbdb="diameter" in SBDB_SAMPLE.get("phys_par", {}),
        notes="NeoWs gives range; SBDB single value",
    ),
    ParameterStatus(
        name="Geometric albedo",
        source="PAIR §2.1",
        available_neows=False,
        available_sbdb="albedo" in SBDB_SAMPLE.get("phys_par", {}),
        notes="SBDB `phys_par.albedo`",
    ),
    ParameterStatus(
        name="Bulk density",
        source="PAIR §2.1",
        available_neows=False,
        available_sbdb="density" in SBDB_SAMPLE.get("phys_par", {}),
        notes="SBDB optional field",
    ),
    ParameterStatus(
        name="Rotation period",
        source="PAIR §2.2",
        available_neows=False,
        available_sbdb="rot_per" in SBDB_SAMPLE.get("phys_par", {}),
        notes="SBDB `phys_par.rot_per`",
    ),
    ParameterStatus(
        name="Spectral class",
        source="PAIR §2.2",
        available_neows=False,
        available_sbdb="spec" in SBDB_SAMPLE.get("phys_par", {}),
        notes="Taxonomy as strength proxy",
    ),
    ParameterStatus(
        name="Relative velocity",
        source="PAIR §2.1",
        available_neows=bool(NEOWS_SAMPLE["close_approach_data"][0]["relative_velocity"]["kilometers_per_second"]),
        available_sbdb=False,
        notes="NeoWs `close_approach_data.relative_velocity`",
    ),
    ParameterStatus(
        name="Miss distance",
        source="PAIR §2.1",
        available_neows=bool(NEOWS_SAMPLE["close_approach_data"][0]["miss_distance"]["kilometers"]),
        available_sbdb=False,
        notes="NeoWs `miss_distance`",
    ),
    ParameterStatus(
        name="Orbital elements (a,e,i,Ω,ω,M)",
        source="PAIR §2.1",
        available_neows=False,
        available_sbdb=set(SBDB_SAMPLE.get("orbit", {})) >= {"a", "e", "i", "om", "w", "ma"},
        notes="Query NeoWs lookup or SBDB",
    ),
]


def maybe_live_fetch(designation: str) -> Optional[Dict[str, object]]:
    if os.getenv("LIVE") != "1":
        return None
    url = "https://ssd-api.jpl.nasa.gov/sbdb.api"
    base_params = {"orb": "1", "phys": "1"}

    def _request(params):
        resp = requests.get(url, params=params, timeout=10)
        resp.raise_for_status()
        return resp.json()

    # Try by designation (DES) first when numeric, otherwise use search string.
    attempts = []
    if designation.isdigit():
        attempts.append({**base_params, "des": designation})
        attempts.append({**base_params, "s": designation})
    else:
        attempts.append({**base_params, "s": designation})
        # If text contains numeric ID, also try DES variant
        tokens = designation.split()
        for token in tokens:
            if token.isdigit():
                attempts.append({**base_params, "des": token})
                break

    last_exc: Optional[Exception] = None
    for params in attempts:
        try:
            return _request(params)
        except requests.RequestException as exc:
            last_exc = exc
            continue

    if last_exc is not None:
        raise last_exc
    raise requests.HTTPError(f"SBDB lookup failed for designation '{designation}'")


def format_bool(value: bool) -> str:
    return "✔" if value else "✘"


def main() -> None:
    print("PAIR parameter availability (NeoWs vs SBDB)")
    print("")
    header = f"{'Parameter':35} {'Source':11} {'NeoWs':6} {'SBDB':6} Notes"
    print(header)
    print("-" * len(header))
    for item in PARAMETERS:
        print(
            f"{item.name:35} {item.source:11} {format_bool(item.available_neows):6} {format_bool(item.available_sbdb):6} {item.notes}"
        )

    target_designation = os.getenv("SBDB_ID", "99942 Apophis")
    live_data = None
    if os.getenv("LIVE") == "1":
        try:
            live_data = maybe_live_fetch(target_designation)
        except requests.RequestException as exc:
            print(f"\nSBDB request failed: {exc}")

    if live_data:
        print("\nLive SBDB fetch:")
        print(json.dumps({
            "designation": live_data.get("object", {}).get("fullname"),
            "albedo": live_data.get("phys_par", {}).get("albedo"),
            "density": live_data.get("phys_par", {}).get("density"),
            "rot_per": live_data.get("phys_par", {}).get("rot_per"),
            "spec": live_data.get("phys_par", {}).get("spec"),
        }, indent=2))
    else:
        print("\nSet LIVE=1 to attempt a real SBDB fetch (network permitting).")


if __name__ == "__main__":
    main()
