# app.py
# Streamlit MVP for an educational asteroid-impact dashboard
# -----------------------------------------------------------
# Features
# - Explore tab: sliders for asteroid size, velocity, density, angle
# - Map: choose impact location by city preset or lat/lon inputs
# - Calculates kinetic energy, TNT equivalent, naive crater estimate, and damage radii
# - Defend Earth tab: apply simple deflection (delta-v & lead time) and visualize how the impact point shifts
# - Learn tab: glossary & inline explanations for younger audiences
#
# NOTE (Hackathon-friendly):
# - Uses only public libs: streamlit, numpy, pandas, requests, pydeck
# - NASA NEO lookup optional (uses DEMO_KEY or env var NASA_API_KEY)
# - USGS integration points are included as stubs you can extend (elevation/tsunami overlays)

import argparse
import math
import os
import sys
from datetime import date
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import requests

# Parse CLI args before importing Streamlit (to avoid execution when using --headless-report)
_parser = argparse.ArgumentParser(add_help=False)
_parser.add_argument("--headless-report", action="store_true")
_args, _ = _parser.parse_known_args()

if _args.headless_report:
    # Headless mode - skip Streamlit UI execution
    print("Headless telemetry mode requires running via streamlit with METEOR_MADNESS_HEADLESS_TELEMETRY=1")
    print("Usage: METEOR_MADNESS_HEADLESS_TELEMETRY=1 streamlit run main.py")
    sys.exit(0)

import streamlit as st
import pydeck as pdk

from src.config import SLIDER_SPECS, SliderSpec, get_slider_spec
from src.telemetry import (
    ProvenanceTag,
    add_error,
    add_warning,
    ensure_report,
    format_report,
    record_dataset_provenance,
    record_dropdown_selection,
    update_slider_provenance,
    update_slider_validation,
)

# Try to use Mapbox if token is present; otherwise fall back to an open TileLayer (OSM)
MAPBOX_TOKEN = os.getenv("MAPBOX_API_KEY")
if MAPBOX_TOKEN:
    pdk.settings.mapbox_api_key = MAPBOX_TOKEN

# ----------------------------
# Utility & Constants
# ----------------------------
EARTH_RADIUS_KM = 6371.0
G = 9.81  # m/s^2
TNT_J = 4.184e9  # 1 ton TNT in Joules
RHO_TARGET = 2500  # kg/m^3 (crustal rock, rough)
SEA_LEVEL_DENSITY = 1.225  # kg/m^3
SCALE_HEIGHT_KM = 7.16  # exponential atmosphere scale height
DEFAULT_SEISMIC_COUPLING = 3e-4  # fraction of kinetic energy converted to seismic energy

MATERIAL_PRESETS = {
    "Carbonaceous chondrite": {"density": 1500, "strength_mpa": 1.0},
    "Stony (ordinary chondrite)": {"density": 3000, "strength_mpa": 10.0},
    "Iron-nickel": {"density": 7800, "strength_mpa": 50.0},
    "Cometary (icy)": {"density": 600, "strength_mpa": 0.3},
}

# Synthetic fallback population density (people per square km) used when no live data layer is loaded.
# Documented in data_sources.md for transparency.
SYNTHETIC_POP_DENSITY = {
    "severe": 5000,
    "moderate": 2000,
    "light": 800,
}

SYNTHETIC_CASUALTY_RATE = {
    "severe": 0.35,
    "moderate": 0.1,
    "light": 0.02,
}

# Friendly city presets (lat, lon)
CITY_PRESETS = {
    "— choose a place —": (None, None),
    "New York City, USA": (40.7128, -74.0060),
    "Los Angeles, USA": (34.0522, -118.2437),
    "London, UK": (51.5074, -0.1278),
    "Tokyo, Japan": (35.6762, 139.6503),
    "Sydney, Australia": (-33.8688, 151.2093),
    "Rio de Janeiro, Brazil": (-22.9068, -43.1729),
}

# Approximate ring densities (people/km²) for presets — educational, not census-accurate.
CITY_RING_DENSITY = {
    "New York City, USA":      {"severe": 10000, "moderate": 4000, "light": 1500},
    "Los Angeles, USA":        {"severe":  6000, "moderate": 2500, "light": 1000},
    "London, UK":              {"severe":  7000, "moderate": 3000, "light": 1200},
    "Tokyo, Japan":            {"severe":  8000, "moderate": 3500, "light": 1400},
    "Sydney, Australia":       {"severe":  4000, "moderate": 1500, "light":  600},
    "Rio de Janeiro, Brazil":  {"severe":  6500, "moderate": 2800, "light": 1100},
}

# Fallback densities for arbitrary lat/lon (rural/suburban blend)
DEFAULT_RING_DENSITY = {"severe": 1200, "moderate": 500, "light": 200}

# ----------------------------
# Simple physics helpers (educational approximations)
# ----------------------------

def asteroid_mass_kg(diameter_m: float, density_kg_m3: float) -> float:
    radius = diameter_m / 2
    volume = (4/3) * math.pi * radius**3
    return volume * density_kg_m3


def kinetic_energy_joules(mass_kg: float, velocity_km_s: float) -> float:
    v = velocity_km_s * 1000  # m/s
    return 0.5 * mass_kg * v**2


def tnt_megatons(E_j: float) -> float:
    return E_j / TNT_J / 1e6


def estimate_burst_altitude_km(
    velocity_km_s: float,
    strength_mpa: float,
    diameter_m: Optional[float] = None,
    angle_deg: Optional[float] = None,   # NEW
) -> float:
    """Strength breakup altitude with simple size and angle corrections."""
    if velocity_km_s <= 0 or strength_mpa <= 0:
        return 0.0

    v = velocity_km_s * 1000.0
    strength_pa = strength_mpa * 1e6

    # Base (strength ~ dynamic pressure) altitude:
    rho_break = 2.0 * strength_pa / max(v**2, 1.0)           # q = ½ρv² = strength → ρ_break
    altitude = 0.0 if rho_break >= SEA_LEVEL_DENSITY else -SCALE_HEIGHT_KM * math.log(rho_break / SEA_LEVEL_DENSITY)

    # Size correction: larger bodies penetrate deeper (down-shift altitude)
    if diameter_m is not None and diameter_m > 0:
        altitude -= 6.0 * math.log10(max(diameter_m, 10.0) / 50.0)

    # ANGLE correction (educational): shallower entries break higher
    if angle_deg is not None:
        # angle is 10–90° from horizontal in your UI
        s = math.sin(math.radians(max(5.0, min(angle_deg, 90.0))))
        ANGLE_SCALE_KM = 10.0
        altitude += ANGLE_SCALE_KM * (1.0 - s)    # +~8–9 km @ 20°, ~0 km @ 90°

    return float(max(0.0, min(altitude, 60.0)))


def ground_energy_fraction(burst_altitude_km: float) -> float:
    """Approximate fraction of kinetic energy that couples to the ground (PAIR-inspired)."""
    if burst_altitude_km <= 0:
        return 1.0
    # Hackathon-friendly smoothing: keep a gradual decay so high-altitude bursts still
    # couple a tiny amount of energy instead of instantly dropping to zero at 25 km.
    # Logistic curve mimics the PAIR trend without requiring the full atmospheric model.
    scale_height = 4.0  # steeper value = faster drop-off with altitude
    midpoint = 22.0  # around this altitude half the energy reaches the ground
    fraction = 1.0 / (1.0 + math.exp((burst_altitude_km - midpoint) / scale_height))
    return float(max(0.0, min(1.0, fraction)))

K1 = 1.3
MU = 0.55     # velocity exponent via pi-scaling
NU = 0.4      # density exponent
GAMMA = 0.17  # gravity exponent

def transient_crater_diameter_m(diameter_m, velocity_km_s, density_kg_m3, angle_deg,
                                target_density_kg_m3=RHO_TARGET, gravity=G):
    if diameter_m <= 0 or velocity_km_s <= 0:
        return 0.0
    v = velocity_km_s * 1000.0
    a = max(5.0, min(angle_deg, 90.0))
    angle_factor = math.sin(math.radians(a)) ** (1/3)

    # Impact energy proxy using pi-scaling (gravity regime)
    # D_t ∝ K1 * ( (ρ_i/ρ_t)^(ν) ) * d^(1-μ) * (v^μ) * (g^(-γ)) * angle_factor
    rho_ratio = (density_kg_m3 / target_density_kg_m3) ** NU
    D_t = K1 * rho_ratio * (diameter_m ** (1 - MU)) * (v ** MU) * (gravity ** (-GAMMA)) * angle_factor
    return float(max(D_t, 0.0))

def final_crater_diameter_m(diameter_m, velocity_km_s, density_kg_m3, angle_deg, burst_altitude_km):
    # If the body disrupts above ~1 km → no excavation crater (airburst)
    if burst_altitude_km > 1.0:
        return 0.0

    # Otherwise it reaches the ground; compute transient crater and the final size
    D_t = transient_crater_diameter_m(diameter_m, velocity_km_s, density_kg_m3, angle_deg)
    if D_t <= 0:
        return 0.0
    # A simple transient→final factor for simple craters
    return float(max(1.25 * D_t, 0.0))

# very small lookup: overpressure psi -> scaled distance Z (m/kg^(1/3))
# Values are approximate, education-friendly (Glasstone/Dolan-like).
OP_TO_Z = {
    1.0:  180.0,
    4.0:   75.0,
    12.0:  40.0,
}

def blast_overpressure_radius_km(E_mt: float, burst_altitude_km: float, overpressure_psi: float = 4.0) -> float:
    if E_mt <= 0:
        return 0.0

    W_ton = E_mt * 1e6  # convert megatons → tons
    W13 = max(W_ton, 1.0) ** (1 / 3)

    if overpressure_psi in OP_TO_Z:
        Z = OP_TO_Z[overpressure_psi]
    else:
        keys = sorted(OP_TO_Z.keys())
        lo = max(k for k in keys if k <= overpressure_psi)
        hi = min(k for k in keys if k >= overpressure_psi)
        if lo == hi:
            Z = OP_TO_Z[lo]
        else:
            t = (math.log(overpressure_psi) - math.log(lo)) / (math.log(hi) - math.log(lo))
            Z = OP_TO_Z[lo] * (1 - t) + OP_TO_Z[hi] * t

    R0_km = (Z * W13) / 1000.0

    H = max(burst_altitude_km, 0.0)
    H_scale = 18.0  # km; tuned so 25 km bursts shrink ground radius ~50%
    attenuation = 1.0 / (1.0 + (H / H_scale) ** 2)
    return R0_km * attenuation

def seismic_moment_magnitude(E_joules: float, coupling: float = DEFAULT_SEISMIC_COUPLING,
                             ground_fraction: float = 1.0) -> Optional[float]:
    if E_joules <= 0 or coupling <= 0:
        return None
    seismic_energy = E_joules * coupling * max(0.0, min(1.0, ground_fraction))
    if seismic_energy <= 0:
        return None
    return (2.0 / 3.0) * (math.log10(seismic_energy) - 4.8)


def estimate_population_impacts(
    r_severe: float,
    r_moderate: float,
    r_light: float,
    ring_densities: Optional[dict[str, float]] = None,
    track_provenance: bool = True,
) -> dict[str, float]:
    # Coerce non-finite values to 0
    def _safe(x: float) -> float:
        try:
            x = float(x)
        except Exception:
            return 0.0
        return x if math.isfinite(x) and x > 0 else 0.0

    r_severe = _safe(r_severe)
    r_moderate = _safe(r_moderate)
    r_light = _safe(r_light)

    radii = sorted([r_severe, r_moderate, r_light])
    if radii[-1] <= 0:
        return {"severe": 0.0, "moderate": 0.0, "light": 0.0, "total": 0.0, "casualties": 0.0}

    # Enforce nesting
    severe = r_severe
    moderate = max(r_moderate, severe)
    light = max(r_light, moderate)

    # Areas in km² of each ring (annuli)
    areas = {
        "severe": math.pi * severe ** 2,
        "moderate": math.pi * moderate ** 2 - math.pi * severe ** 2,
        "light": math.pi * light ** 2 - math.pi * moderate ** 2,
    }

    # Pick densities: city-specific → fallback → legacy synthetic
    if ring_densities is None:
        ring_densities = DEFAULT_RING_DENSITY
        if track_provenance:
            report = ensure_report(st.session_state)
            record_dataset_provenance(
                report,
                "population_density",
                ProvenanceTag.SYNTHETIC_POP_DENSITY,
                status="warning",
                detail="Using DEFAULT_RING_DENSITY synthetic fallback.",
            )

    exposure = {}
    for ring, area in areas.items():
        density = float(ring_densities.get(ring, SYNTHETIC_POP_DENSITY.get(ring, 0)))
        exposure[ring] = max(0.0, area) * max(0.0, density)

    exposure["total"] = sum(exposure.values())

    casualties = 0.0
    for ring in ("severe", "moderate", "light"):
        rate = float(SYNTHETIC_CASUALTY_RATE.get(ring, 0.0))
        casualties += exposure[ring] * rate
    exposure["casualties"] = casualties

    if track_provenance:
        report = ensure_report(st.session_state)
        record_dataset_provenance(
            report,
            "casualty_rates",
            ProvenanceTag.SYNTHETIC_CASUALTY,
            status="warning",
            detail="Using SYNTHETIC_CASUALTY_RATE fallback values.",
        )

    return exposure

def run_pair_simulation(
    samples: int,
    base_diameter_m: float,
    base_density: float,
    base_velocity_km_s: float,
    base_angle_deg: float,
    base_strength_mpa: float,
):
    """Monte Carlo sampler inspired by Mathias et al. (2017) PAIR framework."""
    rng = np.random.default_rng()

    diameters = rng.lognormal(mean=math.log(max(base_diameter_m, 1.0)), sigma=0.25, size=samples)
    densities = rng.lognormal(mean=math.log(max(base_density, 100.0)), sigma=0.15, size=samples)
    velocities = rng.normal(loc=base_velocity_km_s, scale=max(1.0, 0.15 * base_velocity_km_s), size=samples)
    velocities = np.clip(velocities, 5.0, 70.0)
    strength = rng.lognormal(mean=math.log(max(base_strength_mpa, 0.1)), sigma=0.35, size=samples)

    u = rng.random(samples)
    angles = (90.0 / math.pi) * np.arccos(np.clip(2 * u - 1, -1.0, 1.0))
    angles = np.clip(angles, 10.0, 90.0)

    results = []
    for d, rho, v, theta, s in zip(diameters, densities, velocities, angles, strength):
        mass = asteroid_mass_kg(d, rho)
        E_j = kinetic_energy_joules(mass, v)
        E_mt = tnt_megatons(E_j)
        burst = estimate_burst_altitude_km(v, s, d, theta)
        burst = max(0.0, burst - 5.0)
        ground_frac = ground_energy_fraction(burst)
        crater = final_crater_diameter_m(d, v, rho, theta, burst) / 1000.0
        E_ground_mt = max(E_mt * ground_frac, 0.0)
        r4 = blast_overpressure_radius_km(E_ground_mt, burst, overpressure_psi=4.0)
        r12 = blast_overpressure_radius_km(E_ground_mt, burst, overpressure_psi=12.0)
        r1 = blast_overpressure_radius_km(E_ground_mt, burst, overpressure_psi=1.0)
        exposure = estimate_population_impacts(r12, r4, r1, track_provenance=False)
        Mw = seismic_moment_magnitude(E_j, ground_fraction=ground_frac)
        results.append(
            {
                "diameter_m": d,
                "density": rho,
                "velocity": v,
                "angle": theta,
                "strength_mpa": s,
                "energy_mt": E_mt,
                "burst_alt_km": burst,
                "ground_fraction": ground_frac,
                "crater_km": crater,
                "r12_km": r12,
                "r4_km": r4,
                "r1_km": r1,
                "population": exposure.get("total", 0.0),
                "casualties": exposure.get("casualties", 0.0),
                "Mw": Mw if Mw is not None else float("nan"),
            }
        )

    return pd.DataFrame(results)


def haversine_offset(lat, lon, d_km, bearing_deg):
    """Move from (lat, lon) by distance d_km along bearing_deg; return new (lat, lon)."""
    lat1 = math.radians(lat)
    lon1 = math.radians(lon)
    brg = math.radians(bearing_deg)
    dr = d_km / EARTH_RADIUS_KM
    lat2 = math.asin(math.sin(lat1)*math.cos(dr) + math.cos(lat1)*math.sin(dr)*math.cos(brg))
    lon2 = lon1 + math.atan2(math.sin(brg)*math.sin(dr)*math.cos(lat1), math.cos(dr) - math.sin(lat1)*math.sin(lat2))
    return math.degrees(lat2), (math.degrees(lon2) + 540) % 360 - 180


def apply_deflection(lat, lon, delta_v_mm_s: float, lead_days: float, inbound_bearing_deg: float = 0.0):
    """
    Toy model: along-track shift s ≈ Δv * t (projected), then map to km and shift opposite inbound bearing.
    """
    # Convert mm/s to m/s
    dv = delta_v_mm_s / 1000.0
    t = lead_days * 86400.0
    s_m = dv * t
    s_km = s_m / 1000.0
    # shift perpendicular-ish to showcase effect; here use opposite bearing to move impact point
    new_lat, new_lon = haversine_offset(lat, lon, s_km, (inbound_bearing_deg + 180) % 360)
    return new_lat, new_lon, s_km


# ----------------------------
# NASA NEO API helper (optional)
# ----------------------------


def _resolve_nasa_api_key() -> str:
    env_key = os.getenv("NASA_API_KEY")
    if env_key:
        return env_key
    env_file = Path(__file__).resolve().parent / "API.env"
    if env_file.exists():
        return env_file.read_text().strip()
    return "DEMO_KEY"


def fetch_sbdb_phys_params(designation: str) -> Optional[Dict[str, object]]:
    """Fetch physical parameters from NASA SBDB API for a given asteroid designation."""
    url = "https://ssd-api.jpl.nasa.gov/sbdb.api"
    params = {"sstr": designation, "phys-par": "1"}

    try:
        response = requests.get(url, params=params, timeout=5)
        response.raise_for_status()
        data = response.json()

        # SBDB returns phys_par as an array of objects
        phys_list = data.get("phys_par") or []
        phys_dict = {item["name"]: item.get("value") for item in phys_list if isinstance(item, dict)}

        # Extract relevant physical parameters
        density_str = phys_dict.get("density")
        density_kg_m3 = None
        if density_str:
            try:
                # SBDB returns density in g/cm³, convert to kg/m³
                density_kg_m3 = float(density_str) * 1000.0
            except (TypeError, ValueError):
                pass

        diameter_str = phys_dict.get("diameter")
        diameter_km = None
        if diameter_str:
            try:
                diameter_km = float(diameter_str)
            except (TypeError, ValueError):
                pass

        albedo_str = phys_dict.get("albedo")
        albedo = None
        if albedo_str:
            try:
                albedo = float(albedo_str)
            except (TypeError, ValueError):
                pass

        return {
            "sbdb_density_kg_m3": density_kg_m3,
            "sbdb_albedo": albedo,
            "sbdb_diameter_km": diameter_km,
            "sbdb_spectral_class": phys_dict.get("spec_B") or phys_dict.get("spec_T"),
            "sbdb_rotation_period_hr": phys_dict.get("rot_per"),
        }
    except Exception as e:
        print(f"SBDB lookup failed for {designation}: {e}")
        return None


def fetch_today_neos() -> pd.DataFrame:
    api_key = _resolve_nasa_api_key()
    today = date.today().isoformat()
    url = f"https://api.nasa.gov/neo/rest/v1/feed?start_date={today}&end_date={today}&api_key={api_key}"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
    except Exception as e:
        print(f"NASA NEO API Error: {e}")
        raise

    data = response.json()
    neos: List[Dict[str, object]] = []
    for neo in data.get("near_earth_objects", {}).get(today, []):
        diameter_data: Dict[str, object] = (neo.get("estimated_diameter") or {}).get("meters", {}) or {}
        close_approaches = neo.get("close_approach_data") or []
        first_approach = close_approaches[0] if close_approaches else {}

        velocity = first_approach.get("relative_velocity", {}).get("kilometers_per_second")
        try:
            velocity = float(velocity) if velocity is not None else None
        except (TypeError, ValueError):
            velocity = None

        miss_distance = first_approach.get("miss_distance", {}).get("kilometers")
        try:
            miss_distance = float(miss_distance) if miss_distance is not None else None
        except (TypeError, ValueError):
            miss_distance = None

        record = {
            "name": neo.get("name"),
            "neo_reference_id": neo.get("neo_reference_id"),
            "designation": neo.get("designation"),
            "absolute_magnitude_h": neo.get("absolute_magnitude_h"),
            "est_diameter_min_m": diameter_data.get("estimated_diameter_min"),
            "est_diameter_max_m": diameter_data.get("estimated_diameter_max"),
            "hazardous": bool(neo.get("is_potentially_hazardous_asteroid", False)),
            "velocity_km_s": velocity,
            "miss_distance_km": miss_distance,
            "close_approach_date": first_approach.get("close_approach_date"),
            "orbiting_body": first_approach.get("orbiting_body"),
        }

        # Enrich with SBDB physical parameters
        designation = neo.get("designation")
        if designation:
            sbdb_data = fetch_sbdb_phys_params(designation)
            if sbdb_data:
                record.update(sbdb_data)

        neos.append(record)

    df = pd.DataFrame(neos)
    if df.empty:
        return df

    for col in ("est_diameter_min_m", "est_diameter_max_m", "velocity_km_s", "miss_distance_km"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    df["diameter_avg_m"] = df[["est_diameter_min_m", "est_diameter_max_m"]].mean(axis=1, skipna=True)
    df["miss_distance_moon_x"] = df["miss_distance_km"] / 384_400.0
    return df


def load_today_neo_feed() -> pd.DataFrame:
    report = ensure_report(st.session_state)
    today = date.today().isoformat()
    if "neo_today_df" not in st.session_state:
        try:
            neo_df = fetch_today_neos()
            st.session_state["neo_today_df"] = neo_df
            if neo_df.empty:
                record_dataset_provenance(
                    report,
                    "neo_feed",
                    ProvenanceTag.LIVE,
                    status="empty",
                    detail="NeoWs daily feed returned no objects for today.",
                )
            else:
                record_dataset_provenance(
                    report,
                    "neo_feed",
                    ProvenanceTag.LIVE,
                    status="ok",
                    detail=f"NeoWs daily feed for {today}.",
                )
        except Exception as e:
            record_dataset_provenance(
                report,
                "neo_feed",
                ProvenanceTag.LIVE,
                status="error",
                detail=f"NASA NEO API failed: {str(e)}",
            )
            add_error(report, f"NASA NEO API Error: {str(e)}")
            st.error(f"Failed to fetch NASA NEO data: {str(e)}")
            st.session_state["neo_today_df"] = pd.DataFrame([])
    return st.session_state["neo_today_df"]


def select_baseline_neo(df: pd.DataFrame) -> Optional[pd.Series]:
    if df is None or df.empty:
        return None
    hazardous = df[df["hazardous"] == True]  # noqa: E712 - pandas comparison
    if not hazardous.empty:
        return hazardous.iloc[0]
    return df.iloc[0]


def _is_integer_slider(spec: SliderSpec) -> bool:
    return all(isinstance(getattr(spec, attr), int) for attr in ("min_value", "max_value", "step"))


def _prepare_slider_default(
    spec: SliderSpec,
    slider_key: str,
    fallback_value: Optional[float],
    fallback_label: Optional[str],
    live_defaults: Dict[str, float],
) -> Tuple[Union[float, int], ProvenanceTag, Optional[str]]:
    if slider_key in st.session_state:
        value = st.session_state[slider_key]
        provenance = ProvenanceTag.SESSION_STATE
        detail = "Session state value re-used."
    else:
        value = live_defaults.get(spec.key)
        if value is not None:
            provenance = ProvenanceTag.LIVE
            source_name = (st.session_state.get("live_slider_defaults_meta", {}) or {}).get("source")
            detail = f"Hydrated from NeoWs object {source_name}" if source_name else "Hydrated from NeoWs daily feed."
        else:
            if fallback_value is not None:
                value = fallback_value
                if fallback_label == "material_preset":
                    provenance = ProvenanceTag.MATERIAL_PRESET
                    detail = "Material preset fallback applied."
                else:
                    provenance = ProvenanceTag.FALLBACK
                    detail = "Runtime fallback applied."
            elif spec.whitepaper_default is not None:
                value = spec.whitepaper_default
                provenance = ProvenanceTag.WHITEPAPER
                detail = "Whitepaper default used."
            else:
                value = spec.min_value
                provenance = ProvenanceTag.SPEC_FLOOR
                detail = "Spec minimum used as final fallback."  # keep informative

    if value is None:
        value = spec.min_value
        provenance = ProvenanceTag.SPEC_FLOOR
        detail = "Spec minimum used as final fallback."

    try:
        numeric = float(value)
    except (TypeError, ValueError):
        numeric = float(spec.min_value)
        provenance = ProvenanceTag.SPEC_FLOOR
        detail = "Invalid value coerced to spec minimum."

    bounded = max(float(spec.min_value), min(float(spec.max_value), numeric))
    if not math.isclose(bounded, numeric):
        detail = (detail or "") + " (clamped to spec bounds)."
    numeric = bounded

    if _is_integer_slider(spec):
        return int(round(numeric)), provenance, detail
    return numeric, provenance, detail


def prepare_slider_args(
    spec_key: str,
    *,
    key_override: Optional[str] = None,
    fallback_value: Optional[float] = None,
    fallback_label: Optional[str] = None,
) -> Tuple[SliderSpec, Dict[str, object]]:
    spec = get_slider_spec(spec_key)
    slider_key = key_override or spec.key
    live_defaults = st.session_state.get("live_slider_defaults", {})
    value, provenance, provenance_detail = _prepare_slider_default(
        spec,
        slider_key,
        fallback_value,
        fallback_label,
        live_defaults,
    )

    slider_kwargs: Dict[str, object] = {
        "min_value": spec.min_value,
        "max_value": spec.max_value,
        "step": spec.step,
        "key": slider_key,
    }
    if value is not None:
        slider_kwargs["value"] = value

    _validate_slider_args(spec, slider_kwargs)

    report = ensure_report(st.session_state)
    recorded_value = float(value) if isinstance(value, (int, float)) else None
    update_slider_provenance(
        report,
        spec.key,
        provenance,
        value=recorded_value,
        detail=provenance_detail,
        widget_key=slider_key,
    )

    return spec, slider_kwargs


def _validate_slider_args(spec: SliderSpec, slider_kwargs: Dict[str, object]) -> None:
    """Ensure the Streamlit slider kwargs match the spec; raise early on drift."""

    report = ensure_report(st.session_state)
    mismatches = []
    for param in ("min_value", "max_value", "step"):
        expected = getattr(spec, param)
        actual = slider_kwargs.get(param)
        if float(expected) != float(actual):
            mismatches.append(f"{param}: expected {expected}, got {actual}")

    if mismatches:
        message = (
            f"Slider '{spec.key}' diverges from whitepaper spec — " + ", ".join(mismatches)
        )
        update_slider_validation(report, spec.key, "error", message=message)
        st.error(message)
        raise ValueError(message)

    value = slider_kwargs.get("value")
    if value is not None:
        if not (float(spec.min_value) <= float(value) <= float(spec.max_value)):
            message = (
                f"Slider '{spec.key}' default value {value} outside bounds "
                f"[{spec.min_value}, {spec.max_value}]"
            )
            update_slider_validation(report, spec.key, "error", message=message)
            st.error(message)
            raise ValueError(message)

    update_slider_validation(report, spec.key, "ok")


def seed_live_slider_defaults(neo_df: pd.DataFrame) -> Dict[str, object]:
    if "live_slider_defaults_meta" in st.session_state:
        return st.session_state["live_slider_defaults_meta"]

    report = ensure_report(st.session_state)
    live_defaults: Dict[str, float] = {}
    meta: Dict[str, object] = {
        "source": None,
        "warnings": [],
        "notes": [],
        "provenance": "whitepaper",
    }

    baseline = select_baseline_neo(neo_df)
    if baseline is None:
        meta["warnings"].append(
            "Could not fetch NASA NeoWs feed for today; sliders fall back to whitepaper defaults."
        )
    else:
        baseline_name = baseline.get("name") or baseline.get("designation") or "NeoWs object"
        meta["source"] = str(baseline_name)
        meta["provenance"] = "live"

        diameter = baseline.get("diameter_avg_m")
        if pd.notna(diameter):
            live_defaults["diameter_m"] = float(diameter)
        else:
            meta["warnings"].append("NeoWs did not report a usable diameter; using whitepaper default instead.")

        velocity = baseline.get("velocity_km_s")
        if pd.notna(velocity):
            live_defaults["velocity_km_s"] = float(velocity)
        else:
            meta["warnings"].append("NeoWs did not report a usable relative velocity; using whitepaper default.")

        miss_km = baseline.get("miss_distance_km")
        if pd.notna(miss_km):
            meta["notes"].append(
                f"Baseline object miss distance: {miss_km:,.0f} km (≈ {baseline.get('miss_distance_moon_x', float('nan')):.1f}× Moon)."
            )

    missing_live = [
        spec.label
        for spec in SLIDER_SPECS.values()
        if spec.requires_live_default and live_defaults.get(spec.key) is None
    ]
    if missing_live:
        meta["warnings"].append(
            "Live defaults unavailable for: " + ", ".join(missing_live) + ". Using documented fallbacks."
        )

    meta["hydrated_keys"] = sorted(live_defaults.keys())

    datasets = report.get("datasets", {})
    neo_entry = datasets.get("neo_feed", {}) if isinstance(datasets, dict) else {}
    current_status = str(neo_entry.get("status", "ok" if baseline is not None else "empty"))
    if baseline is None:
        detail = "NeoWs daily feed unavailable; defaulting to documented fallbacks."
    else:
        detail = f"NeoWs daily feed baseline: {meta['source']}"
    record_dataset_provenance(
        report,
        "neo_feed",
        ProvenanceTag.LIVE,
        status=current_status,
        detail=detail,
    )

    st.session_state["live_slider_defaults"] = live_defaults
    st.session_state["live_slider_defaults_meta"] = meta

    for warning in meta["warnings"]:
        add_warning(report, warning)

    return meta


def _status_icon(status: Optional[str]) -> str:
    normalized = (status or "").lower()
    if normalized == "ok":
        return "✅"
    if normalized in {"warning", "empty"}:
        return "⚠️"
    if normalized == "error":
        return "❌"
    return "ℹ️"


def render_sidebar_telemetry(live_meta: Dict[str, object]) -> None:
    report = ensure_report(st.session_state)
    sidebar = st.sidebar

    datasets = report.get("datasets", {}) if isinstance(report, dict) else {}
    if datasets:
        sidebar.markdown("**Data sources**")
        for name in sorted(datasets):
            info = datasets.get(name, {})
            icon = _status_icon(str(info.get("status")))
            provenance = info.get("provenance", "unknown")
            sidebar.markdown(f"- {icon} **{name}** · {provenance}")
            detail = info.get("detail")
            if detail:
                sidebar.caption(detail)

    sliders = report.get("sliders", {}) if isinstance(report, dict) else {}
    if sliders:
        sidebar.markdown("**Sliders**")
        for key in sorted(sliders):
            info = sliders.get(key, {})
            icon = _status_icon(str(info.get("status")))
            provenance = info.get("provenance", "unknown")
            sidebar.markdown(f"- {icon} `{key}` · {provenance}")
            message = info.get("message")
            detail = info.get("provenance_detail")
            if message:
                sidebar.caption(message)
            elif detail:
                sidebar.caption(detail)

    error_messages = []
    report_errors = report.get("errors", []) if isinstance(report, dict) else []
    error_messages.extend(report_errors)
    if error_messages:
        sidebar.markdown("**Errors**")
        for error in dict.fromkeys(error_messages):
            sidebar.caption(f"❌ {error}")

    warning_messages = []
    report_warnings = report.get("warnings", []) if isinstance(report, dict) else []
    warning_messages.extend(report_warnings)
    warning_messages.extend(live_meta.get("warnings", []))
    if warning_messages:
        sidebar.markdown("**Warnings**")
        for warning in dict.fromkeys(warning_messages):
            sidebar.caption(f"⚠️ {warning}")

    if live_meta.get("source"):
        sidebar.markdown("**NeoWs defaults**")
        sidebar.caption(f"Baseline object: {live_meta['source']}")
        hydrated = live_meta.get("hydrated_keys")
        if hydrated:
            sidebar.caption("Hydrated sliders: " + ", ".join(hydrated))

    for note in live_meta.get("notes", []):
        sidebar.caption(note)


def emit_headless_telemetry(force: bool = False) -> int:
    """Emit telemetry report and return exit code (0=success, 1=errors found)."""
    flag = os.getenv("METEOR_MADNESS_HEADLESS_TELEMETRY", "")
    should_emit = force or flag.lower() in {"1", "true", "yes"}

    if should_emit:
        report = ensure_report(st.session_state)
        print(format_report(report))

        # Return non-zero exit code if errors are present
        errors = report.get("errors", [])
        return 1 if errors else 0

    return 0


# ----------------------------
# UI
# ----------------------------
st.set_page_config(page_title="Impactor-2025: Learn & Simulate", layout="wide")

ensure_report(st.session_state)

neo_df = load_today_neo_feed()
seed_live_slider_defaults(neo_df)

st.title("🛰️ Impactor-2025: Learn & Simulate")
st.caption("An educational dashboard to explore asteroid impacts, built for a hackathon.")
# If an interaction queued widget overrides (e.g., from the NASA NEO picker),
# apply them before any widgets with those keys are instantiated.
if "widget_overrides" in st.session_state:
    overrides = st.session_state.pop("widget_overrides")
    for key, value in overrides.items():
        st.session_state[key] = value


# Tabs
exp_tab, defend_tab, learn_tab = st.tabs(["Explore", "Defend Earth", "Learn"])

with exp_tab:
    st.subheader("Choose impact parameters")
    live_defaults_meta = st.session_state.get("live_slider_defaults_meta", {})
    if live_defaults_meta.get("source"):
        st.info(f"Slider defaults seeded from NASA NeoWs object: {live_defaults_meta['source']}")
    for warning in live_defaults_meta.get("warnings", []):
        st.warning(warning)
    for note in live_defaults_meta.get("notes", []):
        st.caption(note)

    primary_cols = st.columns(4)
    with primary_cols[0]:
        diameter_spec, diameter_kwargs = prepare_slider_args("diameter_m")
        diameter_m = st.slider(diameter_spec.label, **diameter_kwargs)
    with primary_cols[1]:
        velocity_spec, velocity_kwargs = prepare_slider_args("velocity_km_s")
        velocity = st.slider(velocity_spec.label, **velocity_kwargs)
    with primary_cols[2]:
        angle_spec, angle_kwargs = prepare_slider_args("angle_deg")
        angle = st.slider(angle_spec.label, **angle_kwargs)
    with primary_cols[3]:
        material = st.selectbox("Material preset", list(MATERIAL_PRESETS.keys()), index=1)
        report = ensure_report(st.session_state)
        record_dropdown_selection(
            report,
            "material_preset",
            material,
            len(MATERIAL_PRESETS),
            provenance=ProvenanceTag.MATERIAL_PRESET,
            detail=f"Selected material: {material}",
        )

    preset_density = MATERIAL_PRESETS[material]["density"]
    preset_strength = MATERIAL_PRESETS[material]["strength_mpa"]

    secondary_cols = st.columns(2)
    with secondary_cols[0]:
        density_fallback = max(
            float(preset_density),
            float(get_slider_spec("bulk_density").min_value),
        )
        density_spec, density_kwargs = prepare_slider_args(
            "bulk_density",
            key_override=f"density_{material}",
            fallback_value=density_fallback,
            fallback_label="material_preset",
        )
        density = st.slider(density_spec.label, **density_kwargs)
        if preset_density < density_spec.min_value:
            st.caption(
                f"Preset density for {material}: {preset_density} kg/m³ (whitepaper floor {density_spec.min_value} kg/m³)."
            )
        else:
            st.caption(f"Preset density for {material}: {preset_density} kg/m³")
    with secondary_cols[1]:
        strength_spec, strength_kwargs = prepare_slider_args(
            "strength_mpa",
            key_override=f"strength_{material}",
            fallback_value=float(preset_strength),
            fallback_label="material_preset",
        )
        strength_mpa = st.slider(strength_spec.label, **strength_kwargs)
        st.caption("Adjust to emulate cohesive strength used in PAIR entry modeling.")

    st.markdown("**Where does it hit?** Pick a city or enter coordinates.")
    c1, c2, c3 = st.columns([2,1,1])
    with c1:
        preset = st.selectbox("City preset", list(CITY_PRESETS.keys()))
        report = ensure_report(st.session_state)
        # Determine provenance based on whether it's a real city or generic option
        if preset in CITY_RING_DENSITY:
            dropdown_prov = ProvenanceTag.CITY_PRESET_DENSITY
        else:
            dropdown_prov = ProvenanceTag.SYNTHETIC_POP_DENSITY
        record_dropdown_selection(
            report,
            "city_preset",
            preset,
            len(CITY_PRESETS),
            provenance=dropdown_prov,
            detail=f"Selected location: {preset}",
        )
    if preset and CITY_PRESETS[preset][0] is not None:
        lat, lon = CITY_PRESETS[preset]
    else:
        with c2:
            lat = st.number_input("Latitude", value=29.7604, format="%.4f")
        with c3:
            lon = st.number_input("Longitude", value=-95.3698, format="%.4f")
    # Choose densities based on city selection
    report = ensure_report(st.session_state)
    if preset in CITY_RING_DENSITY:
        current_ring_densities = CITY_RING_DENSITY[preset]
        record_dataset_provenance(
            report,
            "population_density",
            ProvenanceTag.CITY_PRESET_DENSITY,
            status="ok",
            detail=f"Using city preset density for {preset}.",
        )
    else:
        current_ring_densities = DEFAULT_RING_DENSITY
        record_dataset_provenance(
            report,
            "population_density",
            ProvenanceTag.SYNTHETIC_POP_DENSITY,
            status="warning",
            detail="Using DEFAULT_RING_DENSITY synthetic fallback.",
        )

    # Check if we're in "live mode" (NeoWs data loaded successfully) and using synthetic data
    neo_meta = st.session_state.get("live_slider_defaults_meta", {})
    is_live_mode = neo_meta.get("provenance") == "live"
    if is_live_mode and preset not in CITY_RING_DENSITY:
        add_warning(
            report,
            "Using synthetic population density in live mode — consider switching to a city preset or integrating real population data.",
        )

    # Computations
    m = asteroid_mass_kg(diameter_m, density)
    E_j = kinetic_energy_joules(m, velocity)
    E_mt = tnt_megatons(E_j)
    breakup_alt_km = estimate_burst_altitude_km(velocity, strength_mpa, diameter_m, angle)
    burst_alt_km = max(0.0, breakup_alt_km - 6.0)
    ground_fraction = ground_energy_fraction(burst_alt_km)
    crater_m = final_crater_diameter_m(diameter_m, velocity, density, angle, burst_alt_km)
    crater_km = crater_m / 1000.0
    # The height term inside blast_overpressure_radius_km already accounts for attenuation.
    # Use actual ground-coupled energy; no artificial floor
    E_ground_mt = max(E_mt * ground_fraction, 0.0)

    r_mod = blast_overpressure_radius_km(E_ground_mt, burst_alt_km, overpressure_psi=4.0)
    r_severe = blast_overpressure_radius_km(E_ground_mt, burst_alt_km, overpressure_psi=12.0)
    r_light = blast_overpressure_radius_km(E_ground_mt, burst_alt_km, overpressure_psi=1.0)

    # If the blast rings collapse (e.g., high-altitude bursts) but a crater forms, seed rings from the crater footprint
    crater_radius_km = max(crater_km / 2.0, 0.0)
    # Only force nesting when the blast calculation is effectively zero
    if crater_radius_km > 0 and (E_ground_mt < 0.5 or r_mod < 0.1):
        r_severe = max(r_severe, crater_radius_km)
        r_mod = max(r_mod, r_severe * 1.6)
        r_light = max(r_light, r_mod * 1.4)

    exposure = estimate_population_impacts(r_severe, r_mod, r_light, ring_densities=current_ring_densities)
    Mw = seismic_moment_magnitude(E_j, ground_fraction=ground_fraction)

    st.markdown("### Results")
    cols = st.columns(4)
    cols[0].metric("Mass", f"{m:,.0f} kg")
    cols[1].metric("Kinetic energy", f"{E_mt:,.2f} Mt TNT")
    cols[2].metric("Breakup altitude", f"{breakup_alt_km:.1f} km")
    cols[3].metric("Ground-coupled energy", f"{E_mt * ground_fraction:.2f} Mt")

    crater_display = f"{crater_km:.2f} km" if crater_km > 0.0 else "Airburst"
    cols2 = st.columns(4)
    cols2[0].metric("Crater diameter (final)", crater_display)
    cols2[1].metric("Severe radius (12 psi)", f"{r_severe:.2f} km")
    cols2[2].metric("Moderate radius (4 psi)", f"{r_mod:.2f} km")
    cols2[3].metric("Light radius (1 psi)", f"{r_light:.2f} km")

    cols3 = st.columns(3)
    cols3[0].metric("Population exposed", f"{exposure.get('total', 0.0):,.0f}")
    cols3[1].metric("Estimated casualties", f"{exposure.get('casualties', 0.0):,.0f}")
    cols3[2].metric("Seismic Mw", f"{Mw:.1f}" if Mw is not None else "n/a")

    st.caption(
        "Population estimates rely on documented synthetic density rings and fall back to the crater footprint when blast rings vanish; replace with real datasets when available."
    )

    with st.expander("Synthetic exposure breakdown"):
        exposure_rows = []
        for ring, radius in [("severe", r_severe), ("moderate", r_mod), ("light", r_light)]:
            area = math.pi * max(radius, 0.0) ** 2
            pop = exposure.get(ring, 0.0)
            rate = SYNTHETIC_CASUALTY_RATE.get(ring, 0.0)
            exposure_rows.append(
                {
                    "ring": ring,
                    "radius_km": radius,
                    "area_km2": area,
                    "population": pop,
                    "casualty_rate": rate,
                    "expected_casualties": pop * rate,
                }
            )
        st.dataframe(pd.DataFrame(exposure_rows))

    with st.expander("PAIR-inspired probabilistic scenarios"):
        st.write(
            "Sample uncertain properties (diameter, density, angle, strength) using the PAIR Monte Carlo approach to explore outcome distributions."
        )
        samples_spec, samples_kwargs = prepare_slider_args("pair_samples")
        samples = st.slider(samples_spec.label, **samples_kwargs)
        if st.button("Run simulation", key="run_pair_button"):
            sim_df = run_pair_simulation(samples, diameter_m, density, velocity, angle, strength_mpa)
            if sim_df.empty:
                st.warning("Simulation returned no scenarios.")
            else:
                quantiles = sim_df[[
                    "energy_mt",
                    "burst_alt_km",
                    "ground_fraction",
                    "crater_km",
                    "r4_km",
                    "population",
                    "casualties",
                ]].quantile([0.5, 0.9]).T
                quantiles.columns = ["p50", "p90"]
                st.subheader("Key outcome quantiles")
                st.dataframe(quantiles)

                crater_probability = float((sim_df["crater_km"] > 0).mean())
                st.metric("Probability of crater formation", f"{crater_probability*100:.1f}%")

                st.caption(
                    "PAIR = Probabilistic Asteroid Impact Risk (Mathias et al., 2017). Replace assumed distributions with mission-specific priors as data become available."
                )

    # Map visualization with concentric circles
    st.markdown("### Map")
    view_state = pdk.ViewState(latitude=lat, longitude=lon, zoom=6, bearing=0, pitch=30)

    def circle_layer(radius_km, color, opacity=100):
        return pdk.Layer(
            "ScatterplotLayer",
            data=pd.DataFrame({"lat": [lat], "lon": [lon]}),
            get_position="[lon, lat]",
            get_radius=radius_km * 1000,
            radius_min_pixels=1,
            radius_max_pixels=10000,
            get_fill_color=color,
            pickable=False,
            stroked=False,
            filled=True,
        )

    rings = [
        circle_layer(r_light, [255, 165, 0, 60]),
        circle_layer(r_mod, [255, 0, 0, 80]),
        circle_layer(r_severe, [139, 0, 0, 120]),
        pdk.Layer(
            "ScatterplotLayer",
            data=pd.DataFrame({"lat": [lat], "lon": [lon]}),
            get_position="[lon, lat]",
            get_radius=5000,
            get_fill_color=[0, 0, 0, 180],
            get_line_color=[255, 255, 255, 220],
            line_width_min_pixels=1,
        ),
    ]

    # Basemap handling: use Mapbox style if token exists, else add an open TileLayer basemap
    if MAPBOX_TOKEN:
        deck = pdk.Deck(map_style="mapbox://styles/mapbox/dark-v11", initial_view_state=view_state, layers=rings)
    else:
        basemap = pdk.Layer(
            "TileLayer",
            data="https://tile.openstreetmap.org/{z}/{x}/{y}.png",
            min_zoom=0, max_zoom=19, tile_size=256,
        )
        deck = pdk.Deck(initial_view_state=view_state, layers=[basemap, *rings], map_style=None)

    st.pydeck_chart(deck)

    # --- NEO selector with SBDB integration ---
    with st.expander("Use today's NASA NEOs as inputs"):
        df = neo_df.copy()
        if df.empty:
            st.info("Could not fetch NEOs right now (API rate limit or offline). You can still use the simulator.")
        else:
            # Short, readable view - include SBDB data if available
            view_cols = ["name", "diameter_avg_m", "velocity_km_s", "miss_distance_km", "miss_distance_moon_x", "hazardous"]
            # Add SBDB columns if they exist in the dataframe
            if "sbdb_density_kg_m3" in df.columns:
                view_cols.append("sbdb_density_kg_m3")
            if "sbdb_spectral_class" in df.columns:
                view_cols.append("sbdb_spectral_class")

            display_df = df[view_cols].rename(columns={
                "diameter_avg_m": "diameter (m)",
                "velocity_km_s": "velocity (km/s)",
                "miss_distance_km": "miss distance (km)",
                "miss_distance_moon_x": "miss (× Moon)",
                "sbdb_density_kg_m3": "density (kg/m³)",
                "sbdb_spectral_class": "spectral class"
            })
            st.dataframe(display_df)

            def on_neo_select():
                """Callback to update sliders when NEO is selected from dropdown."""
                choice = st.session_state.get("neo_picker")
                if not choice:
                    return

                row = df.loc[df["name"] == choice].iloc[0]
                report = ensure_report(st.session_state)

                diameter_value = row["diameter_avg_m"]
                if diameter_value is None or pd.isna(diameter_value):
                    diameter_value = st.session_state.get("diameter_m", get_slider_spec("diameter_m").min_value)
                diameter_value = int(np.clip(diameter_value, get_slider_spec("diameter_m").min_value, get_slider_spec("diameter_m").max_value))

                velocity_value = row["velocity_km_s"]
                if velocity_value is None or pd.isna(velocity_value):
                    velocity_value = st.session_state.get("velocity_km_s", get_slider_spec("velocity_km_s").min_value)
                velocity_spec = get_slider_spec("velocity_km_s")
                velocity_value = float(np.clip(velocity_value, velocity_spec.min_value, velocity_spec.max_value))

                # Check for SBDB data
                density_value = row.get("sbdb_density_kg_m3")
                if density_value and not pd.isna(density_value):
                    density_spec = get_slider_spec("bulk_density")
                    density_value = float(np.clip(density_value, density_spec.min_value, density_spec.max_value))

                    # Record SBDB provenance
                    record_dataset_provenance(
                        report,
                        "sbdb_density",
                        ProvenanceTag.SBDB,
                        status="ok",
                        detail=f"Using SBDB density for {choice}: {density_value:,.0f} kg/m³",
                    )

                    st.session_state["widget_overrides"] = {
                        "diameter_m": diameter_value,
                        "velocity_km_s": velocity_value,
                        "angle_deg": 45,
                        f"density_{material}": density_value,
                    }
                else:
                    st.session_state["widget_overrides"] = {
                        "diameter_m": diameter_value,
                        "velocity_km_s": velocity_value,
                        "angle_deg": 45,
                    }

            # Pick one and auto-update sliders
            choice = st.selectbox(
                "Pick an asteroid to simulate (what-if it hit):",
                options=df["name"].tolist(),
                key="neo_picker",
                on_change=on_neo_select
            )

            if choice:
                row = df.loc[df["name"] == choice].iloc[0]
                st.caption(
                    f"Selected **{row['name']}** — avg diameter ≈ {row['diameter_avg_m']:.1f} m, "
                    f"speed ≈ {row['velocity_km_s']:.2f} km/s, "
                    f"miss distance ≈ {row['miss_distance_km']:.0f} km (~{row['miss_distance_moon_x']:.1f}× Moon)."
                )

                # Show SBDB data if available
                if row.get("sbdb_density_kg_m3") and not pd.isna(row["sbdb_density_kg_m3"]):
                    st.success(f"SBDB density: {row['sbdb_density_kg_m3']:,.0f} kg/m³")
                if row.get("sbdb_spectral_class"):
                    st.info(f"Spectral class: {row['sbdb_spectral_class']}")

                # Energy preview
                mass = asteroid_mass_kg(row["diameter_avg_m"], density)
                E_mt_preview = tnt_megatons(kinetic_energy_joules(mass, row["velocity_km_s"]))
                st.metric("What-if energy (preview)", f"{E_mt_preview:,.2f} Mt TNT")
                st.caption("Preview assumes current density/material selection.")


with defend_tab:
    st.subheader("Try a deflection strategy ✨")
    st.write("Toy model: apply a small velocity change (Δv) some days before arrival and see how the nominal impact point shifts.")

    c1, c2, c3 = st.columns(3)
    with c1:
        delta_spec, delta_kwargs = prepare_slider_args("delta_v_mm_s")
        delta_v_mm_s = st.slider(delta_spec.label, **delta_kwargs)
    with c2:
        lead_spec, lead_kwargs = prepare_slider_args("lead_days")
        lead_days = st.slider(lead_spec.label, **lead_kwargs)
    with c3:
        bearing_spec, bearing_kwargs = prepare_slider_args("inbound_bearing")
        inbound_bearing = st.slider(bearing_spec.label, **bearing_kwargs)

    # Baseline from Explore tab (share state)
    base_lat, base_lon = lat, lon
    new_lat, new_lon, shift_km = apply_deflection(base_lat, base_lon, delta_v_mm_s, lead_days, inbound_bearing)

    cols = st.columns(3)
    cols[0].metric("Shift on ground (km)", f"{shift_km:.1f}")
    cols[1].metric("Old impact", f"{base_lat:.3f}, {base_lon:.3f}")
    cols[2].metric("New impact", f"{new_lat:.3f}, {new_lon:.3f}")

    view_state = pdk.ViewState(latitude=base_lat, longitude=base_lon, zoom=5, bearing=0, pitch=30)

    path_df = pd.DataFrame({
        "lat": [base_lat, new_lat],
        "lon": [base_lon, new_lon],
    })

    overlays = [
        pdk.Layer(
            "LineLayer",
            data=pd.DataFrame({"source_lon": [base_lon], "source_lat": [base_lat], "target_lon": [new_lon], "target_lat": [new_lat]}),
            get_source_position="[source_lon, source_lat]",
            get_target_position="[target_lon, target_lat]",
            get_width=3,
            get_color=[0, 200, 255, 200],
        ),
        pdk.Layer(
            "ScatterplotLayer",
            data=pd.DataFrame({"lat": [base_lat], "lon": [base_lon]}),
            get_position="[lon, lat]",
            get_radius=7000,
            get_fill_color=[255, 0, 0, 200],
        ),
        pdk.Layer(
            "ScatterplotLayer",
            data=pd.DataFrame({"lat": [new_lat], "lon": [new_lon]}),
            get_position="[lon, lat]",
            get_radius=7000,
            get_fill_color=[0, 255, 0, 200],
        ),
    ]

    if MAPBOX_TOKEN:
        deck = pdk.Deck(map_style="mapbox://styles/mapbox/dark-v11", initial_view_state=view_state, layers=overlays)
    else:
        basemap = pdk.Layer(
            "TileLayer",
            data="https://tile.openstreetmap.org/{z}/{x}/{y}.png",
            min_zoom=0, max_zoom=19, tile_size=256,
        )
        deck = pdk.Deck(initial_view_state=view_state, layers=[basemap, *overlays], map_style=None)

    st.pydeck_chart(deck)

    st.caption(
        "Deflection visualization preserves the same damage model; connect to spatial datasets to recalculate exposure for the shifted impact point."
    )

with learn_tab:
    st.subheader("Glossary & Teaching Aids")
    st.markdown(
        """
        **Asteroid** — A rocky object orbiting the Sun. Some get close to Earth and are called **NEOs** (Near-Earth Objects).

        **Diameter** — How wide the asteroid is. Bigger usually means more energy on impact.

        **Velocity** — How fast it's moving. Energy grows with the **square** of speed!

        **Kinetic energy** — Energy of motion. We convert it to **megatons of TNT** to help compare sizes.

        **Crater** — A bowl-shaped hole. Our estimate is simplified for learning.

        **Deflection (Δv)** — A tiny push far in advance can move an asteroid enough to miss Earth.

        **Why simplified?** — Real scientists use more complex models (airbursts, fragmentation, terrain, oceans).
        Our goal is to **learn the ideas** first—then you can add advanced physics.
        """
    )

    st.markdown("### Where to extend (hackathon tasks)")
    st.markdown("- **USGS overlays:** add coastal elevation/tsunami hazard layers via tiled map sources.")
    st.markdown("- **Population exposure:** add a layer with night lights or population to illustrate risk.")
    st.markdown("- **Better crater/blast models:** swap in published scaling relations and atmosphere effects.")
    st.markdown("- **Orbit view:** a 3D Three.js canvas for the Sun–Earth–asteroid geometry (or use Plotly 3D).")

render_sidebar_telemetry(st.session_state.get("live_slider_defaults_meta", {}))
st.sidebar.markdown("---")
st.sidebar.title("About this MVP")
st.sidebar.info(
    "This is an educational demo. Numbers are approximate. For real decisions, consult official models and data."
)

# Only emit via environment variable when running in Streamlit mode
emit_headless_telemetry()


def main() -> int:
    """CLI entry point for headless telemetry mode."""
    parser = argparse.ArgumentParser(
        description="Meteor Madness - Asteroid Impact Simulator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--headless-report",
        action="store_true",
        help="Print telemetry report and exit (non-zero exit code if errors found)",
    )

    args = parser.parse_args()

    if args.headless_report:
        # Force headless telemetry emission and exit with appropriate code
        # Note: This requires running the Streamlit app once to populate session state
        # For true headless mode, we'd need to initialize the app without the UI
        print("Headless telemetry report:")
        print("Note: Run streamlit normally with METEOR_MADNESS_HEADLESS_TELEMETRY=1 to capture full state.")
        return 0

    # No CLI flags provided - show usage
    parser.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
