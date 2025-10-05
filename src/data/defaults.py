from __future__ import annotations

import os
from typing import Any, Dict, Optional, Tuple

import requests

from src.api.nasa_neo import (
    NeoWsClient,
    NeoWsError,
    estimate_mean_diameter_km,
    extract_orbital_elements,
)


DEFAULT_NEO_ID = os.getenv("NEO_ID", "3542519")  # 99942 Apophis
DEFAULT_SBDB_ID = os.getenv("SBDB_ID", "99942 Apophis")
SBDB_API_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api"


class SBDBError(RuntimeError):
    """Raised when the SBDB API returns an HTTP or parsing error."""


MATERIAL_PRESETS = {
    "carbonaceous_chondrite": {"density": 1500, "strength_mpa": 1.0},
    "stony_chondrite": {"density": 3000, "strength_mpa": 10.0},
    "iron_nickel": {"density": 7800, "strength_mpa": 50.0},
    "cometary_icy": {"density": 600, "strength_mpa": 0.3},
}


TAXONOMY_TO_PRESET = {
    "A": "stony_chondrite",
    "B": "carbonaceous_chondrite",
    "C": "carbonaceous_chondrite",
    "D": "cometary_icy",
    "E": "iron_nickel",
    "F": "carbonaceous_chondrite",
    "G": "carbonaceous_chondrite",
    "I": "iron_nickel",
    "K": "stony_chondrite",
    "L": "stony_chondrite",
    "M": "iron_nickel",
    "O": "stony_chondrite",
    "P": "cometary_icy",
    "Q": "stony_chondrite",
    "R": "stony_chondrite",
    "S": "stony_chondrite",
    "T": "carbonaceous_chondrite",
    "U": "stony_chondrite",
    "V": "stony_chondrite",
    "X": "iron_nickel",
    "Z": "carbonaceous_chondrite",
}


def _safe_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _select_close_approach(neo_payload: Dict[str, Any]) -> Tuple[Optional[float], Optional[float], Optional[str], Optional[Dict[str, Any]]]:
    """Return preferred close-approach velocity/miss-distance plus provenance."""

    close_approaches = (neo_payload.get("close_approach_data") or [])
    if not close_approaches:
        return None, None, None, None

    selected = None
    source_label: Optional[str] = None
    for approach in close_approaches:
        if approach.get("orbiting_body") == "Earth":
            selected = approach
            source_label = "NeoWs close_approach_data (Earth)"
            break

    if selected is None:
        selected = close_approaches[0]
        source_label = "NeoWs close_approach_data (first entry)"

    rel_vel = _safe_float((selected.get("relative_velocity") or {}).get("kilometers_per_second"))
    miss_km = _safe_float((selected.get("miss_distance") or {}).get("kilometers"))

    return rel_vel, miss_km, source_label, selected


def material_from_taxonomy(spec: Optional[str]) -> str:
    if not spec:
        return "stony_chondrite"
    key = spec.strip().upper()
    preset = TAXONOMY_TO_PRESET.get(key)
    if not preset and key:
        preset = TAXONOMY_TO_PRESET.get(key[0])
    return preset or "stony_chondrite"


def normalize_sbdb_payload(payload: Dict[str, Any]) -> Dict[str, Any]:
    phys = payload.get("phys_par")
    if isinstance(phys, list):
        flattened: Dict[str, Any] = {}
        for entry in phys:
            if not isinstance(entry, dict):
                continue
            key = entry.get("name") or entry.get("field")
            if not key:
                continue
            flattened[key] = entry.get("value")
        payload["phys_par"] = flattened
    return payload


def fetch_sbdb_payload(designation: str) -> Dict[str, Any]:
    params = {
        "sstr": designation,
        "phys-par": "1",
        "full-prec": "1",
    }
    try:
        response = requests.get(SBDB_API_URL, params=params, timeout=20)
        response.raise_for_status()
    except requests.RequestException as exc:
        raise SBDBError(f"SBDB request failed: {exc}") from exc
    try:
        data = response.json()
    except ValueError as exc:
        raise SBDBError("SBDB response was not valid JSON") from exc
    return normalize_sbdb_payload(data)


def extract_sbdb_phys(payload: Dict[str, Any]) -> Dict[str, Optional[float]]:
    phys = payload.get("phys_par") or {}
    density_raw = _safe_float(phys.get("density"))
    density_kg_m3: Optional[float] = None
    if density_raw is not None:
        density_kg_m3 = density_raw * 1000.0 if density_raw < 100 else density_raw
    rot_period = _safe_float(phys.get("rot_per"))
    albedo = _safe_float(phys.get("albedo"))
    return {
        "density_kg_m3": density_kg_m3,
        "rotation_period_hr": rot_period,
        "albedo": albedo,
        "spectral_class": (phys.get("spec") or "").strip() or None,
    }


def fetch_neows_object(neo_id: str) -> Dict[str, Any]:
    with NeoWsClient() as client:
        return client.get_neo(neo_id)


def resolve_density_strength(
    sbdb_density: Optional[float],
    spectral_class: Optional[str],
) -> Tuple[Optional[float], Optional[float], str]:
    if sbdb_density is not None:
        closest_material = min(
            MATERIAL_PRESETS.items(),
            key=lambda item: abs(item[1]["density"] - sbdb_density),
        )[0]
        return sbdb_density, MATERIAL_PRESETS[closest_material]["strength_mpa"], "SBDB density"

    material = material_from_taxonomy(spectral_class)
    preset = MATERIAL_PRESETS[material]
    return preset["density"], preset["strength_mpa"], f"Taxonomy {spectral_class or 'unknown'}"


def get_reference_defaults(
    neo_id: str = DEFAULT_NEO_ID,
    sbdb_id: str = DEFAULT_SBDB_ID,
) -> Dict[str, Any]:
    defaults: Dict[str, Any] = {
        "diameter_m": 150,
        "velocity_km_s": 18.0,
        "material": "stony_chondrite",
        "density": MATERIAL_PRESETS["stony_chondrite"]["density"],
        "strength_mpa": MATERIAL_PRESETS["stony_chondrite"]["strength_mpa"],
        "source": "synthetic fallback",
    }

    field_sources: Dict[str, str] = {
        "diameter_m": "Synthetic fallback (150 m default)",
        "velocity_km_s": "Synthetic fallback (18 km/s default)",
        "density": "Material preset: Stony (ordinary chondrite)",
        "strength_mpa": "Material preset: Stony (ordinary chondrite)",
        "material": "Material preset: Stony (ordinary chondrite)",
        "absolute_magnitude_h": "Not provided (synthetic fallback)",
        "albedo": "Not provided (synthetic fallback)",
        "rotation_period_hr": "Not provided (synthetic fallback)",
        "semi_major_axis_au": "Synthetic placeholder",
        "eccentricity": "Synthetic placeholder",
        "inclination_deg": "Synthetic placeholder",
        "ascending_node_longitude_deg": "Synthetic placeholder",
        "argument_of_periapsis_deg": "Synthetic placeholder",
        "miss_distance_km": "Not provided (synthetic fallback)",
        "taxonomy": "Not provided (synthetic fallback)",
    }

    try:
        neo_payload = fetch_neows_object(neo_id)
    except NeoWsError as exc:
        defaults["source_error"] = str(exc)
        neo_payload = {}

    if neo_payload:
        defaults["neo_name"] = neo_payload.get("name")
        defaults["neo_designation"] = neo_payload.get("designation")
        H_mag = _safe_float(neo_payload.get("absolute_magnitude_h"))
        if H_mag is not None:
            defaults["absolute_magnitude_h"] = H_mag
            field_sources["absolute_magnitude_h"] = "NeoWs absolute_magnitude_h"
        diameter_km = estimate_mean_diameter_km(neo_payload)
        if diameter_km is not None:
            defaults["diameter_m"] = float(diameter_km * 1000.0)
            field_sources["diameter_m"] = "NeoWs estimated_diameter (avg min/max)"
        velocity, miss_distance_km, approach_source, approach_snapshot = _select_close_approach(neo_payload)
        if velocity is not None:
            defaults["velocity_km_s"] = velocity
            if approach_source:
                field_sources["velocity_km_s"] = approach_source
        if miss_distance_km is not None:
            defaults["miss_distance_km"] = miss_distance_km
            field_sources["miss_distance_km"] = "NeoWs close_approach_data miss_distance"
        if approach_snapshot is not None:
            defaults["close_approach_snapshot"] = approach_snapshot
            defaults["close_approach_source_label"] = approach_source

        orbital = extract_orbital_elements(neo_payload)
        if orbital:
            if orbital.get("semi_major_axis_au") is not None:
                defaults["semi_major_axis_au"] = orbital["semi_major_axis_au"]
                field_sources["semi_major_axis_au"] = "NeoWs orbital_data.a"
            if orbital.get("eccentricity") is not None:
                defaults["eccentricity"] = orbital["eccentricity"]
                field_sources["eccentricity"] = "NeoWs orbital_data.e"
            if orbital.get("inclination_deg") is not None:
                defaults["inclination_deg"] = orbital["inclination_deg"]
                field_sources["inclination_deg"] = "NeoWs orbital_data.i"
            if orbital.get("ascending_node_longitude_deg") is not None:
                defaults["ascending_node_longitude_deg"] = orbital["ascending_node_longitude_deg"]
                field_sources["ascending_node_longitude_deg"] = "NeoWs orbital_data.Ω"
            if orbital.get("argument_of_periapsis_deg") is not None:
                defaults["argument_of_periapsis_deg"] = orbital["argument_of_periapsis_deg"]
                field_sources["argument_of_periapsis_deg"] = "NeoWs orbital_data.ω"

    try:
        sbdb_payload = fetch_sbdb_payload(sbdb_id)
    except SBDBError as exc:
        defaults.setdefault("source_error", str(exc))
        sbdb_payload = {}

    if sbdb_payload:
        defaults["sbdb_fullname"] = (
            sbdb_payload.get("object", {}).get("fullname")
            or sbdb_payload.get("object", {}).get("des")
            or sbdb_id
        )
        phys = extract_sbdb_phys(sbdb_payload)
        density, strength, provenance = resolve_density_strength(
            phys.get("density_kg_m3"),
            phys.get("spectral_class"),
        )
        if density is not None:
            defaults["density"] = density
            field_sources["density"] = f"{provenance}"
        if strength is not None:
            defaults["strength_mpa"] = strength
            field_sources["strength_mpa"] = f"{provenance}"
        defaults["material"] = material_from_taxonomy(phys.get("spectral_class"))
        field_sources["material"] = (
            f"SBDB taxonomy {phys.get('spectral_class') or 'unknown'}"
            if phys.get("spectral_class")
            else "Material preset fallback"
        )
        defaults["taxonomy"] = phys.get("spectral_class")
        if phys.get("spectral_class"):
            field_sources["taxonomy"] = "SBDB phys_par.spec"
        defaults["provenance"] = provenance
        defaults["sbdb_density"] = phys.get("density_kg_m3")
        if phys.get("albedo") is not None:
            defaults["albedo"] = phys["albedo"]
            field_sources["albedo"] = "SBDB phys_par.albedo"
        if phys.get("rotation_period_hr") is not None:
            defaults["rotation_period_hr"] = phys["rotation_period_hr"]
            field_sources["rotation_period_hr"] = "SBDB phys_par.rot_per"

    defaults["field_sources"] = field_sources

    return defaults


__all__ = [
    "DEFAULT_NEO_ID",
    "DEFAULT_SBDB_ID",
    "SBDB_API_URL",
    "SBDBError",
    "MATERIAL_PRESETS",
    "TAXONOMY_TO_PRESET",
    "_safe_float",
    "material_from_taxonomy",
    "normalize_sbdb_payload",
    "fetch_sbdb_payload",
    "extract_sbdb_phys",
    "fetch_neows_object",
    "resolve_density_strength",
    "get_reference_defaults",
]
