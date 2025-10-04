"""Cross-check white-paper parameters against live NeoWs and SBDB APIs.

Usage (from repo root):

    python3 tests/run_data_source_check.py

Environment variables:

- `NASA_API_KEY`: overrides DEMO_KEY for NeoWs requests (recommended).
- `NEO_ID`: choose a specific NeoWs object (default: 3542519, Apophis).
- `SBDB_ID`: choose a specific SBDB designation or search string (default: 99942 Apophis).

The script must reach both APIs; any network or credential failure aborts with a
clear message so the operator knows the data are not current.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from typing import Any, Dict, Iterable, Optional, Tuple

import requests


def load_env_if_present(path: str = ".env") -> None:
    if not os.path.exists(path):
        return
    try:
        with open(path, "r", encoding="utf-8") as env_file:
            for raw_line in env_file:
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                if "=" not in line:
                    continue
                key, _, value = line.partition("=")
                key = key.strip()
                value = value.strip().strip('"').strip("'")
                os.environ.setdefault(key, value)
    except OSError:
        pass


load_env_if_present()


@dataclass
class ParameterSpec:
    label: str
    neo_path: Optional[Iterable[Any]]
    sbdb_path: Optional[Iterable[Any]]
    notes: str = ""


PARAMETERS = [
    ParameterSpec("Absolute magnitude (H)", ("absolute_magnitude_h",), ("phys_par", "H")),
    ParameterSpec("Diameter max (m)", ("estimated_diameter", "meters", "estimated_diameter_max"), ("phys_par", "diameter"), "SBDB value in km"),
    ParameterSpec("Geometric albedo", None, ("phys_par", "albedo")),
    ParameterSpec("Bulk density (kg/m^3)", None, ("phys_par", "density")),
    ParameterSpec("Rotation period (hr)", None, ("phys_par", "rot_per")),
    ParameterSpec("Spectral class", None, ("phys_par", "spec")),
    ParameterSpec("Relative velocity (km/s)", ("close_approach_data", 0, "relative_velocity", "kilometers_per_second"), None),
    ParameterSpec("Miss distance (km)", ("close_approach_data", 0, "miss_distance", "kilometers"), None),
    ParameterSpec("Semi-major axis a (AU)", ("orbital_data", "a"), ("orbit", "a")),
    ParameterSpec("Eccentricity e", ("orbital_data", "e"), ("orbit", "e")),
    ParameterSpec("Inclination i (deg)", ("orbital_data", "i"), ("orbit", "i")),
    ParameterSpec("Argument of periapsis ω (deg)", ("orbital_data", "omega"), ("orbit", "w")),
    ParameterSpec("Longitude ascending node Ω (deg)", ("orbital_data", "ascending_node_longitude"), ("orbit", "om")),
]


# Allowed SBDB lookup params from https://ssd-api.jpl.nasa.gov/doc/sbdb.html.
# Filtering avoids viewer-only flags (e.g., orb/log) that trigger parameter errors.
ALLOWED_SBDB_PARAMS = {
    "sstr",
    "des",
    "spk",
    "full-prec",
    "phys-par",
    "alt-des",
    "alt-spk",
    "alt-orbits",
    "cov",
    "ca-data",
    "ca-time",
    "ca-tunc",
    "nv-fmt",
    "discovery",
    "orbit-defs",
    "r-notes",
    "r-observer",
    "cd-epoch",
    "cd-tp",
}


def normalize_sbdb_payload(data: Dict[str, Any]) -> Dict[str, Any]:
    """Flatten sections that come back as name/value lists for easier lookup."""

    orbit = data.get("orbit")
    if isinstance(orbit, dict):
        elements = orbit.get("elements")
        if isinstance(elements, list):
            for element in elements:
                if not isinstance(element, dict):
                    continue
                name = element.get("name")
                value = element.get("value")
                if name and name not in orbit:
                    orbit[name] = value

    phys_par = data.get("phys_par")
    if isinstance(phys_par, list):
        flattened: Dict[str, Any] = {}
        for entry in phys_par:
            if not isinstance(entry, dict):
                continue
            key = entry.get("name") or entry.get("field")
            if not key:
                continue
            flattened[key] = entry.get("value")
        data["phys_par"] = flattened

    return data


def get_nested(data: Dict[str, Any], path: Optional[Iterable[Any]]) -> Optional[Any]:
    if path is None:
        return None
    curr: Any = data
    for key in path:
        try:
            curr = curr[key]
        except (KeyError, IndexError, TypeError):
            return None
    return curr


def fetch_neows(neo_id: str) -> Tuple[Dict[str, Any], Optional[str]]:
    api_key = os.getenv("NASA_API_KEY", "DEMO_KEY")
    url = f"https://api.nasa.gov/neo/rest/v1/neo/{neo_id}?api_key={api_key}"
    try:
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
        return resp.json(), None
    except requests.RequestException as exc:
        raise RuntimeError(f"NeoWs request failed: {exc}") from exc


def build_sbdb_attempts(designation: str) -> Iterable[Dict[str, Any]]:
    base = {"full-prec": "1", "phys-par": "1"}
    tokens = designation.split()
    numeric = next((tok for tok in tokens if tok.isdigit()), None)

    def with_base(extra: Dict[str, Any]) -> Dict[str, Any]:
        combined = {**base, **extra}
        return {key: value for key, value in combined.items() if key in ALLOWED_SBDB_PARAMS}

    if numeric:
        yield with_base({"des": numeric})

        try:
            spk_id = str(2_000_000 + int(numeric))
        except ValueError:
            spk_id = None
        if spk_id:
            yield with_base({"spk": spk_id})

    # Always try the full designation first to match exact names like "99942 Apophis".
    yield with_base({"sstr": designation})

    # If we have multiple tokens, try the non-numeric tail (e.g., "Apophis").
    if tokens:
        tail_tokens = [tok for tok in tokens if not tok.isdigit()]
        for tail in tail_tokens:
            yield with_base({"sstr": tail})


def fetch_sbdb(designation: str) -> Tuple[Dict[str, Any], Optional[str]]:
    url = "https://ssd-api.jpl.nasa.gov/sbdb.api"
    errors: list[str] = []
    for params in build_sbdb_attempts(designation):
        try:
            resp = requests.get(url, params=params, timeout=15)
            resp.raise_for_status()
            return normalize_sbdb_payload(resp.json()), None
        except requests.HTTPError as exc:
            body = None
            if exc.response is not None:
                try:
                    body = exc.response.text
                except Exception:  # pragma: no cover - best effort debug aid
                    body = None
            if body:
                errors.append(f"{exc}; response body: {body}")
            else:
                errors.append(str(exc))
            continue
        except requests.RequestException as exc:
            errors.append(str(exc))
            continue

    raise RuntimeError(f"SBDB request failed after attempts: {'; '.join(errors)}")


def format_cell(value: Optional[Any]) -> str:
    if value is None:
        return "✘"
    if isinstance(value, str):
        text = value
    else:
        text = str(value)
    if len(text) > 10:
        text = text[:9] + "…"
    return f"✔ {text}"


def main() -> None:
    neo_id = os.getenv("NEO_ID", "3542519")
    sbdb_id = os.getenv("SBDB_ID", "99942 Apophis")

    try:
        neo_data, neo_note = fetch_neows(neo_id)
        sbdb_data, sbdb_note = fetch_sbdb(sbdb_id)
    except RuntimeError as exc:
        print(f"ERROR: {exc}")
        sys.exit(2)

    print(f"NeoWs object: {neo_data.get('name', neo_id)}")
    print(f"SBDB object: {sbdb_data.get('object', {}).get('fullname', sbdb_id)}")
    print("")

    header = f"{'Parameter':30} {'NeoWs':15} {'SBDB':15} Notes"
    print(header)
    print("-" * len(header))

    for spec in PARAMETERS:
        neo_value = get_nested(neo_data, spec.neo_path)
        sbdb_value = get_nested(sbdb_data, spec.sbdb_path)
        print(
            f"{spec.label:30} {format_cell(neo_value):15} {format_cell(sbdb_value):15} {spec.notes}"
        )

    footer_msgs = [note for note in (neo_note, sbdb_note) if note]
    if footer_msgs:
        print("\nWarnings:")
        for msg in footer_msgs:
            print(f"- {msg}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
