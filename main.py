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

import math
import os
from datetime import date
from typing import Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd
import streamlit as st
import pydeck as pdk

from src.i18n import t, set_lang, get_lang, AVAILABLE_LANGS

from src.api.nasa_neo import NeoWsClient, NeoWsError
from src.data.defaults import (
    DEFAULT_NEO_ID,
    DEFAULT_SBDB_ID,
    MATERIAL_PRESETS,
    SBDBError,
    _safe_float,
    extract_sbdb_phys,
    fetch_sbdb_payload,
    get_reference_defaults,
    material_from_taxonomy,
    resolve_density_strength,
    fetch_neows_object,
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
# External data helpers
# ----------------------------


@st.cache_data(show_spinner=False)
def load_reference_defaults() -> Dict[str, Any]:
    return get_reference_defaults()


@st.cache_data(show_spinner=False)
def cached_sbdb_payload(designation: str) -> Dict[str, Any]:
    return fetch_sbdb_payload(designation)


@st.cache_data(show_spinner=False)
def cached_neo_payload(neo_id: str) -> Dict[str, Any]:
    return fetch_neows_object(neo_id)

_DEFAULTS_TRACE_EMITTED = False

def print_startup_data_trace(
    defaults: Dict[str, Any], 
    control_state: Optional[Dict[str, Any]] = None,
) -> None:
    global _DEFAULTS_TRACE_EMITTED
    if _DEFAULTS_TRACE_EMITTED:
        return
    _DEFAULTS_TRACE_EMITTED = True

    print("\n=== Dashboard data provenance snapshot ===")
    neo_label = defaults.get("neo_name") or defaults.get("neo_designation") or DEFAULT_NEO_ID
    sbdb_label = defaults.get("sbdb_fullname") or DEFAULT_SBDB_ID
    provenance = defaults.get("provenance") or defaults.get("source", "synthetic fallback")
    if defaults.get("source_error"):
        print(f"[warn] data fetch issues: {defaults['source_error']}")
    print(f"NeoWs target : {neo_label} (ID {DEFAULT_NEO_ID})")
    print(f"SBDB target  : {sbdb_label}")
    print(f"Density source: {provenance}")

    field_sources = defaults.get("field_sources") or {}
    if field_sources:
        print("Field-level sources:")
        for key in sorted(field_sources):
            print(f"  - {key}: {field_sources[key]}")

    try:
        from tests.run_data_source_check import PARAMETERS as WHITEPAPER_PARAMS, get_nested
    except ImportError:
        print("[info] Parameter comparison table unavailable (tests package not importable).")
        return

    try:
        neo_payload = cached_neo_payload(DEFAULT_NEO_ID)
        neo_error = None
    except NeoWsError as exc:
        neo_payload = {}
        neo_error = str(exc)
    except Exception as exc:  # pragma: no cover - defensive log
        neo_payload = {}
        neo_error = str(exc)

    try:
        sbdb_payload = cached_sbdb_payload(DEFAULT_SBDB_ID)
        sbdb_error = None
    except SBDBError as exc:
        sbdb_payload = {}
        sbdb_error = str(exc)
    except Exception as exc:  # pragma: no cover - defensive log
        sbdb_payload = {}
        sbdb_error = str(exc)

    if neo_error:
        print(f"[warn] NeoWs payload unavailable: {neo_error}")
    if sbdb_error:
        print(f"[warn] SBDB payload unavailable: {sbdb_error}")

    selected_approach = defaults.get("close_approach_snapshot") or {}
    selected_velocity = defaults.get("velocity_km_s")
    selected_miss = defaults.get("miss_distance_km")

    header = f"{ 'Parameter':35} { 'NeoWs value':>18} { 'SBDB value':>18} { 'Dashboard':>18} { 'Source':>18}"
    print("\n" + header)
    print("-" * len(header))

    dashboard_field_map = {
        "Absolute magnitude (H)": defaults.get("absolute_magnitude_h"),
        "Diameter max (m)": defaults.get("diameter_m"),
        "Bulk density (kg/m^3)": defaults.get("density"),
        "Relative velocity (km/s)": defaults.get("velocity_km_s"),
        "Spectral class": defaults.get("taxonomy"),
        "Geometric albedo": defaults.get("albedo"),
        "Rotation period (hr)": defaults.get("rotation_period_hr"),
        "Miss distance (km)": defaults.get("miss_distance_km"),
        "Semi-major axis a (AU)": defaults.get("semi_major_axis_au"),
        "Eccentricity e": defaults.get("eccentricity"),
        "Inclination i (deg)": defaults.get("inclination_deg"),
        "Argument of periapsis ω (deg)": defaults.get("argument_of_periapsis_deg"),
        "Longitude ascending node Ω (deg)": defaults.get("ascending_node_longitude_deg"),
    }

    source_key_map = {
        "Absolute magnitude (H)": "absolute_magnitude_h",
        "Diameter max (m)": "diameter_m",
        "Bulk density (kg/m^3)": "density",
        "Relative velocity (km/s)": "velocity_km_s",
        "Spectral class": "taxonomy",
        "Geometric albedo": "albedo",
        "Rotation period (hr)": "rotation_period_hr",
        "Miss distance (km)": "miss_distance_km",
        "Semi-major axis a (AU)": "semi_major_axis_au",
        "Eccentricity e": "eccentricity",
        "Inclination i (deg)": "inclination_deg",
        "Argument of periapsis ω (deg)": "argument_of_periapsis_deg",
        "Longitude ascending node Ω (deg)": "ascending_node_longitude_deg",
    }

    def _fmt(value: object) -> str:
        if value is None or value == "":
            return "—"
        if isinstance(value, float):
            return f"{value:,.3g}"
        if isinstance(value, int):
            return f"{value:,.0f}"
        text = str(value)
        return text if len(text) <= 17 else text[:14] + "…"

    for spec in WHITEPAPER_PARAMS:
        if spec.label == "Relative velocity (km/s)" and selected_velocity is not None:
            neo_value = selected_velocity
        elif spec.label == "Miss distance (km)" and selected_miss is not None:
            neo_value = selected_miss
        else:
            neo_value = get_nested(neo_payload, spec.neo_path) if neo_payload else None
        sbdb_value = get_nested(sbdb_payload, spec.sbdb_path) if sbdb_payload else None
        dashboard_value = dashboard_field_map.get(spec.label, "not used")
        source_label = "—"
        source_key = source_key_map.get(spec.label)
        if source_key and source_key in field_sources and dashboard_value != "not used":
            source_label = field_sources[source_key]

        print(
            f"{spec.label:35} "
            f"{_fmt(neo_value):>18} "
            f"{_fmt(sbdb_value):>18} "
            f"{_fmt(dashboard_value):>18} "
            f"{_fmt(source_label):>18}"
        )

    if control_state:
        print("\nUI controls at startup:")
        for label, value in control_state.items():
            print(f"  - {label}: {value}")


def gather_ui_control_state() -> Dict[str, Any]:
    """Collect current slider/selectbox inputs for console transparency."""

    control_map = [
        (t("labels.diameter"), "diameter_m"),
        (t("labels.velocity"), "velocity_km_s"),
        (t("labels.angle"), "angle_deg"),
        (t("labels.material"), "material_preset"),
        (t("labels.density"), "bulk_density"),
        (t("labels.strength"), "bulk_strength"),
        (t("app.city_preset"), "city_preset"),
        (t("snapshot.impact_lat_used"), "impact_lat_used"),
        (t("snapshot.impact_lon_used"), "impact_lon_used"),
        (t("labels.latitude"), "impact_lat_manual"),
        (t("labels.longitude"), "impact_lon_manual"),
        (t("app.samples_label"), "monte_carlo_samples"),
        (t("deflect.delta_v_label"), "deflect_delta_v"),
        (t("deflect.lead_time_label"), "deflect_lead_days"),
        (t("deflect.bearing_label"), "deflect_bearing"),
        (t("app.neo_pick_prompt"), "neo_pick"),
    ]

    state: Dict[str, Any] = {}
    for label, key in control_map:
        if key in st.session_state:
            state[label] = st.session_state.get(key)
    return state

def ensure_session_defaults() -> None:
    if st.session_state.get("_defaults_initialized"):
        return

    defaults = load_reference_defaults()

    diameter_val = defaults.get("diameter_m")
    if diameter_val is not None:
        st.session_state.setdefault("diameter_m", int(round(diameter_val)))
    else:
        st.session_state.setdefault("diameter_m", 150.0)

    velocity_val = defaults.get("velocity_km_s")
    if velocity_val is not None:
        st.session_state.setdefault("velocity_km_s", float(velocity_val))
    else:
        st.session_state.setdefault("velocity_km_s", 18.0)

    st.session_state.setdefault("angle_deg", 45)

    material = defaults.get("material") or "Stony (ordinary chondrite)"
    st.session_state.setdefault("material_preset", material)

    density_val = defaults.get("density")
    if density_val is not None:
        st.session_state.setdefault("bulk_density", float(density_val))
    else:
        st.session_state.setdefault(
            "bulk_density", MATERIAL_PRESETS[material]["density"]
        )

    strength_val = defaults.get("strength_mpa")
    if strength_val is not None:
        st.session_state.setdefault("bulk_strength", float(strength_val))
    else:
        st.session_state.setdefault(
            "bulk_strength", MATERIAL_PRESETS[material]["strength_mpa"]
        )

    st.session_state["defaults_metadata"] = defaults
    st.session_state["_defaults_initialized"] = True
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
        exposure = estimate_population_impacts(r12, r4, r1)
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

@st.cache_data(show_spinner=False)
def fetch_today_neos() -> pd.DataFrame:
    today = date.today().isoformat()
    try:
        with NeoWsClient() as client:
            feed = client.get_feed(today, end_date=today, detailed=False)
    except NeoWsError:
        return pd.DataFrame([])

    neos = []
    objects = (feed.get("near_earth_objects") or {}).get(today, [])
    for neo in objects:
        diameter_info = (neo.get("estimated_diameter") or {}).get("meters") or {}
        close_approach = (neo.get("close_approach_data") or [])
        approach = close_approach[0] if close_approach else {}
        velocity = _safe_float((approach.get("relative_velocity") or {}).get("kilometers_per_second"))
        miss_km = _safe_float((approach.get("miss_distance") or {}).get("kilometers"))
        neos.append(
            {
                "name": neo.get("name"),
                "neo_reference_id": neo.get("neo_reference_id"),
                "designation": neo.get("designation"),
                "absolute_magnitude_h": _safe_float(neo.get("absolute_magnitude_h")),
                "est_diameter_min_m": _safe_float(diameter_info.get("estimated_diameter_min")),
                "est_diameter_max_m": _safe_float(diameter_info.get("estimated_diameter_max")),
                "hazardous": bool(neo.get("is_potentially_hazardous_asteroid", False)),
                "velocity_km_s": velocity,
                "miss_distance_km": miss_km,
            }
        )

    return pd.DataFrame(neos)


# ----------------------------
# UI
# ----------------------------
st.set_page_config(page_title="Impactor-2025: Learn & Simulate", layout="wide")

# Language selector (default to English). Stored in session_state to persist across interactions.
if "lang" not in st.session_state:
    # prefer query param if present
    q = st.experimental_get_query_params().get("lang", [None])[0]
    if q:
        try:
            set_lang(q)
            st.session_state["lang"] = q
        except Exception:
            set_lang("en")
            st.session_state["lang"] = "en"
    else:
        set_lang("en")
        st.session_state["lang"] = "en"

with st.sidebar:
    try:
        langs = AVAILABLE_LANGS or ["en"]
    except Exception:
        langs = ["en"]
    # display local (native) names for languages while storing language codes internally
    # mapping: display name -> code
    lang_display_map = {
        "English": "en",
        "繁體中文": "zh-Hant",
        "Español": "es",
        "Français": "fr",
    }
    # build selectbox options based on available langs, showing native names only
    display_options = [k for k, v in lang_display_map.items() if v in langs]
    # compute current selection display label
    current_display = next((k for k, v in lang_display_map.items() if v == st.session_state["lang"]), "English")
    chosen_display = st.selectbox(t("sidebar.language_label"), display_options, index=display_options.index(current_display) if current_display in display_options else 0)
    chosen = lang_display_map.get(chosen_display, "en")
    if chosen != st.session_state["lang"]:
        set_lang(chosen)
        st.session_state["lang"] = chosen
        # reflect in URL so links share language
        st.experimental_set_query_params(lang=chosen)

st.title(t("app.title"))
st.caption(t("app.caption"))

ensure_session_defaults()
defaults_meta = st.session_state.get("defaults_metadata", {})

# If an interaction queued widget overrides (e.g., from the NASA NEO picker),
# apply them before any widgets with those keys are instantiated.
if "widget_overrides" in st.session_state:
    overrides = st.session_state.pop("widget_overrides")
    for key, value in overrides.items():
        st.session_state[key] = value


# Tabs
exp_tab, defend_tab, learn_tab = st.tabs([t("app.explore") or t("app.explore") or t("explore"), t("app.defend_earth") or t("defend_earth"), t("app.learn_label") or t("learn_label")])

with exp_tab:
    st.subheader(t("app.choose_params"))
    if defaults_meta:
        sbdb_name = defaults_meta.get("sbdb_fullname") or DEFAULT_SBDB_ID
        st.caption(
            t("app.caption")
        )
        field_sources = defaults_meta.get("field_sources") or {}

        def _format_metric(value: Optional[float], digits: int = 3) -> str:
            if value is None:
                return "—"
            if isinstance(value, float):
                return f"{value:,.{digits}g}"
            return str(value)

        metric_specs = [
            (t("labels.abs_magnitude") or "Abs magnitude (H)", "absolute_magnitude_h", 3),
            (t("app.geom_albedo"), "albedo", 3),
            (t("app.rotation_period_hr"), "rotation_period_hr", 3),
            (t("app.miss_distance_km"), "miss_distance_km", 4),
            (t("app.spectral_class"), "taxonomy", 3),
        ]
        metric_cols = st.columns(len(metric_specs))
        for col, (label, key, digits) in zip(metric_cols, metric_specs):
            display_value = defaults_meta.get(key)
            if key == "taxonomy":
                col.metric(label, display_value or "—")
            else:
                col.metric(label, _format_metric(display_value, digits))


        with st.expander(t("expanders.white_paper_inputs"), expanded=False):
                table_specs = [
                    ("Absolute magnitude (H)", "absolute_magnitude_h", "mag"),
                    ("Diameter max (m)", "diameter_m", "m"),
                    (t("app.geom_albedo"), "albedo", ""),
                    ("Bulk density (kg/m^3)", "density", "kg/m³"),
                    (t("app.rotation_period_hr"), "rotation_period_hr", "hr"),
                    (t("app.spectral_class"), "taxonomy", ""),
                    ("Relative velocity (km/s)", "velocity_km_s", "km/s"),
                    (t("app.miss_distance_km"), "miss_distance_km", "km"),
                    ("Semi-major axis a (AU)", "semi_major_axis_au", "AU"),
                    (t("labels.absolute_magnitude") or "Absolute magnitude (H)", "absolute_magnitude_h", "mag"),
                    ("Eccentricity e", "eccentricity", ""),
                    ("Inclination i (deg)", "inclination_deg", "°"),
                    ("Argument of periapsis ω (deg)", "argument_of_periapsis_deg", "°"),
                    ("Longitude ascending node Ω (deg)", "ascending_node_longitude_deg", "°"),
                ]

                def _format_table_value(raw_value: Any) -> str:
                    if raw_value is None or raw_value == "":
                        return "—"
                    if isinstance(raw_value, int):
                        return f"{raw_value:,}"
                    if isinstance(raw_value, float):
                        formatted = f"{raw_value:,.6g}"
                        if formatted.startswith("."):
                            formatted = "0" + formatted
                        elif formatted.startswith("-."):
                            formatted = formatted.replace("-.", "-0.", 1)
                        return formatted
                    return str(raw_value)

                rows = []
                for label, key, unit in table_specs:
                    raw_value = defaults_meta.get(key)
                    rows.append(
                        {
                            "Parameter": label,
                            "Value": _format_table_value(raw_value),
                            "Units": unit,
                            "Source": field_sources.get(key, defaults_meta.get("provenance", "—")),
                        }
                    )

                # use_container_width is the supported way to make the dataframe fill the layout
                st.dataframe(pd.DataFrame(rows), use_container_width=True)
    primary_cols = st.columns(4)
    with primary_cols[0]:
        diameter_default = int(st.session_state.get("diameter_m", 150))
        diameter_kwargs = {
            "min_value": 10,
            "max_value": 2000,
            "step": 10,
            "key": "diameter_m",
        }
        if "diameter_m" not in st.session_state:
            diameter_kwargs["value"] = diameter_default
        diameter_m = st.slider(t("labels.diameter"), **diameter_kwargs)
    with primary_cols[1]:
        velocity_default = float(st.session_state.get("velocity_km_s", 18.0))
        velocity_kwargs = {
            "min_value": 5.0,
            "max_value": 70.0,
            "step": 0.5,
            "key": "velocity_km_s",
        }
        if "velocity_km_s" not in st.session_state:
            velocity_kwargs["value"] = velocity_default
        velocity = st.slider(t("labels.velocity"), **velocity_kwargs)
    with primary_cols[2]:
        angle_default = int(st.session_state.get("angle_deg", 45))
        angle_kwargs = {
            "min_value": 10,
            "max_value": 90,
            "step": 1,
            "key": "angle_deg",
        }
        if "angle_deg" not in st.session_state:
            angle_kwargs["value"] = angle_default
        angle = st.slider(t("labels.angle"), **angle_kwargs)
    with primary_cols[3]:
        material_options = list(MATERIAL_PRESETS.keys())
        default_material = st.session_state.get("material_preset", defaults_meta.get("material", material_options[1]))
        try:
            material_index = material_options.index(default_material)
        except ValueError:
            material_index = 1
        material = st.selectbox(
            t("labels.material"),
            material_options,
            index=material_index,
            key="material_preset",
        )

    preset_density = MATERIAL_PRESETS[material]["density"]
    preset_strength = MATERIAL_PRESETS[material]["strength_mpa"]

    secondary_cols = st.columns(2)
    with secondary_cols[0]:
        density_default = float(st.session_state.get("bulk_density", preset_density))
        density_kwargs = {
            "min_value": 300,
            "max_value": 9000,
            "step": 50,
            "key": "bulk_density",
        }
        if "bulk_density" not in st.session_state:
            density_kwargs["value"] = int(round(density_default))
        density = st.slider(t("labels.density"), **density_kwargs)
        provenance = defaults_meta.get("provenance")
        if provenance:
            st.caption(t("app.default_from_provenance", prov=provenance, dens=f"{density_default:,.0f}"))
        else:
            st.caption(t("app.preset_density_for", material=material, dens=f"{preset_density}"))
    with secondary_cols[1]:
        strength_default = float(st.session_state.get("bulk_strength", preset_strength))
        strength_kwargs = {
            "min_value": 0.1,
            "max_value": 300.0,
            "step": 0.1,
            "key": "bulk_strength",
        }
        if "bulk_strength" not in st.session_state:
            strength_kwargs["value"] = float(strength_default)
        strength_mpa = st.slider(t("labels.strength"), **strength_kwargs)
    st.caption(t("app.adjust_strength_caption"))

    st.markdown("**" + t("app.where_hit") + "**")
    c1, c2, c3 = st.columns([2,1,1])
    with c1:
        preset = st.selectbox(t("app.city_preset"), list(CITY_PRESETS.keys()), key="city_preset")
    if preset and CITY_PRESETS[preset][0] is not None:
        lat, lon = CITY_PRESETS[preset]
    else:
        with c2:
            lat = st.number_input(t("labels.latitude"), value=29.7604, format="%.4f", key="impact_lat_manual")
        with c3:
            lon = st.number_input(t("labels.longitude"), value=-95.3698, format="%.4f", key="impact_lon_manual")
    # Choose densities based on city selection
    if preset in CITY_RING_DENSITY:
        current_ring_densities = CITY_RING_DENSITY[preset]
    else:
        current_ring_densities = DEFAULT_RING_DENSITY

    st.session_state["impact_lat_used"] = float(lat)
    st.session_state["impact_lon_used"] = float(lon)

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

    st.markdown("### " + t("app.results"))
    cols = st.columns(4)
    cols[0].metric(t("metrics.mass"), f"{m:,.0f} kg")
    cols[1].metric(t("metrics.energy"), f"{E_mt:,.2f} Mt TNT")
    cols[2].metric(t("metrics.breakup_alt"), f"{breakup_alt_km:.1f} km")
    cols[3].metric(t("metrics.ground_energy"), f"{E_mt * ground_fraction:.2f} Mt")

    crater_display = f"{crater_km:.2f} km" if crater_km > 0.0 else t("app.airburst")
    cols2 = st.columns(4)
    cols2[0].metric(t("metrics.crater"), crater_display)
    cols2[1].metric(t("metrics.severe_radius"), f"{r_severe:.2f} km")
    cols2[2].metric(t("metrics.moderate_radius"), f"{r_mod:.2f} km")
    cols2[3].metric(t("metrics.light_radius"), f"{r_light:.2f} km")

    cols3 = st.columns(3)
    cols3[0].metric(t("metrics.population_exposed"), f"{exposure.get('total', 0.0):,.0f}")
    cols3[1].metric(t("metrics.estimated_casualties"), f"{exposure.get('casualties', 0.0):,.0f}")
    cols3[2].metric(t("metrics.seismic_mw"), f"{Mw:.1f}" if Mw is not None else t("app.na"))

    st.caption(t("app.population_estimates_note"))

    with st.expander(t("expanders.synthetic_exposure")):
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

    with st.expander(t("expanders.pair_scenarios")):
        st.write(t("app.pair_sample_explain"))
        samples = st.slider(t("app.samples_label"), 100, 1000, 300, step=50, key="monte_carlo_samples")
        if st.button(t("app.run_sim_button"), key="run_pair_button"):
            sim_df = run_pair_simulation(samples, diameter_m, density, velocity, angle, strength_mpa)
            if sim_df.empty:
                st.warning(t("app.simulation_no_scenarios"))
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
                st.subheader(t("app.key_outcome_quantiles") if t("app.key_outcome_quantiles", default=None) else "Key outcome quantiles")
                st.dataframe(quantiles)

                crater_probability = float((sim_df["crater_km"] > 0).mean())
                st.metric(t("metrics.probability_crater"), f"{crater_probability*100:.1f}%")

                st.caption(t("app.pair_note"))

    # Map visualization with concentric circles
    st.markdown("### " + t("app.map"))
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

    # --- PATCH 2: make NEOs usable + educational summaries ---
    with st.expander(t("expanders.use_today_neos")):
        df = fetch_today_neos()
        if df.empty:
            st.info(t("app.neo_fetch_unavailable") or "Could not fetch NEOs right now (API rate limit or offline). You can still use the simulator.")
        else:
            # Derive average diameter and lunar-distance multiples for context
            df = df.copy()
            df["diameter_avg_m"] = 0.5 * (df["est_diameter_min_m"] + df["est_diameter_max_m"])
            MOON_KM = 384_400.0
            df["miss_distance_moon_x"] = df["miss_distance_km"] / MOON_KM

            # Short, readable view
            view_cols = ["name", "diameter_avg_m", "velocity_km_s", "miss_distance_km", "miss_distance_moon_x", "hazardous"]
            st.dataframe(df[view_cols].rename(columns={
                "diameter_avg_m": t("neo.cols.diameter_avg_m") or "diameter_avg_m (m)",
                "miss_distance_km": t("neo.cols.miss_distance_km") or "miss_distance (km)",
                "miss_distance_moon_x": t("neo.cols.miss_distance_moon_x") or "miss distance (× Moon)"
            }))

            # Pick one and push into the simulator
            choice = st.selectbox(
                t("app.neo_pick_prompt"),
                options=df["name"].tolist(),
                key="neo_pick",
            )
            row = df.loc[df["name"] == choice].iloc[0]

            st.caption(
                t(
                    "app.neo_selected_caption",
                    name=row.get("name"),
                    diameter=f"{row['diameter_avg_m']:.1f}",
                    speed=f"{row['velocity_km_s']:.2f}",
                    miss_km=f"{row['miss_distance_km']:.0f}",
                    moonx=f"{row['miss_distance_moon_x']:.1f}",
                )
            )

            colA, colB = st.columns(2)
            with colA:
                if st.button(t("app.use_neo_button") or "Use this NEO in the simulator"):
                    diameter_value = row["diameter_avg_m"]
                    if diameter_value is None or pd.isna(diameter_value):
                        diameter_value = st.session_state["diameter_m"]
                    diameter_value = int(np.clip(diameter_value, 10, 2000))

                    velocity_value = row["velocity_km_s"]
                    if velocity_value is None or pd.isna(velocity_value):
                        velocity_value = st.session_state["velocity_km_s"]
                    velocity_value = float(np.clip(velocity_value, 5.0, 70.0))

                    sbdb_density = None
                    sbdb_strength = None
                    sbdb_material = st.session_state.get("material_preset", "Stony (ordinary chondrite)")
                    sbdb_provenance = None
                    sbdb_taxonomy = None
                    sbdb_fullname = None
                    designation = row.get("designation") or row.get("name")
                    if designation:
                        try:
                            sbdb_payload = cached_sbdb_payload(designation)
                            phys = extract_sbdb_phys(sbdb_payload)
                            sbdb_taxonomy = phys.get("spectral_class")
                            sbdb_material = material_from_taxonomy(sbdb_taxonomy)
                            sbdb_density, sbdb_strength, sbdb_provenance = resolve_density_strength(
                                phys.get("density_kg_m3"),
                                sbdb_taxonomy,
                            )
                            sbdb_fullname = (
                                sbdb_payload.get("object", {}).get("fullname")
                                or designation
                            )
                        except SBDBError as exc:
                            st.warning(t("app.sbdb_lookup_failed", err=str(exc)))

                    if sbdb_density is None:
                        sbdb_density = st.session_state.get("bulk_density")
                    if sbdb_strength is None:
                        sbdb_strength = st.session_state.get("bulk_strength")

                    override_values = {
                        "diameter_m": diameter_value,
                        "velocity_km_s": velocity_value,
                        # Reset to a representative entry angle for clarity when swapping asteroids
                        "angle_deg": 45,
                        "bulk_density": float(sbdb_density) if sbdb_density is not None else st.session_state.get("bulk_density", float(MATERIAL_PRESETS[sbdb_material]["density"])),
                        "bulk_strength": float(sbdb_strength) if sbdb_strength is not None else st.session_state.get("bulk_strength", float(MATERIAL_PRESETS[sbdb_material]["strength_mpa"])),
                        "material_preset": sbdb_material,
                    }
                    st.session_state["widget_overrides"] = override_values
                    st.session_state["defaults_metadata"] = {
                        "neo_name": row.get("name"),
                        "neo_designation": row.get("designation"),
                        "sbdb_fullname": sbdb_fullname or row.get("name"),
                        "material": sbdb_material,
                        "density": override_values["bulk_density"],
                        "strength_mpa": override_values["bulk_strength"],
                        "provenance": sbdb_provenance or "SBDB/NeoWs selection",
                        "taxonomy": sbdb_taxonomy,
                    }
                    st.rerun()
            with colB:
                # A quick educational “scale” card
                mass = asteroid_mass_kg(row["diameter_avg_m"], density)
                E_mt_preview = tnt_megatons(kinetic_energy_joules(mass, row["velocity_km_s"]))
                st.metric(t("metrics.whatif_preview") or "What-if energy (preview)", f"{E_mt_preview:,.2f} Mt TNT")
                st.caption(t("app.preview_caption") or "Preview assumes current density/material selection.")


with defend_tab:
    st.subheader(t("app.deflection_title") or "Try a deflection strategy ✨")
    st.write(t("learn.learn_text") or "Toy model: apply a small velocity change (Δv) some days before arrival and see how the nominal impact point shifts.")

    c1, c2, c3 = st.columns(3)
    with c1:
        delta_v_mm_s = st.slider(t("deflect.delta_v_label") or "Δv (mm/s)", 0.0, 5.0, 0.5, step=0.1, key="deflect_delta_v")
    with c2:
        lead_days = st.slider(t("deflect.lead_time_label") or "Lead time (days)", 0, 3650, 365, step=30, key="deflect_lead_days")
    with c3:
        inbound_bearing = st.slider(t("deflect.bearing_label") or "Inbound bearing (°)", 0, 359, 90, key="deflect_bearing")

    # Baseline from Explore tab (share state)
    base_lat, base_lon = lat, lon
    new_lat, new_lon, shift_km = apply_deflection(base_lat, base_lon, delta_v_mm_s, lead_days, inbound_bearing)

    cols = st.columns(3)
    cols[0].metric(t("metrics.shift_km", default="Shift on ground (km)"), f"{shift_km:.1f}")
    cols[1].metric(t("metrics.old_impact", default="Old impact"), f"{base_lat:.3f}, {base_lon:.3f}")
    cols[2].metric(t("metrics.new_impact", default="New impact"), f"{new_lat:.3f}, {new_lon:.3f}")

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

    st.caption(t("app.deflect_caption") or "Deflection visualization preserves the same damage model; connect to spatial datasets to recalculate exposure for the shifted impact point.")

with learn_tab:
    st.subheader(t("learn.glossary_title") or "Glossary & Teaching Aids")
    # long learn text moved to the i18n catalog
    st.markdown(t("learn.learn_text") or """
        **Asteroid** — A rocky object orbiting the Sun. Some get close to Earth and are called **NEOs** (Near-Earth Objects).

        **Diameter** — How wide the asteroid is. Bigger usually means more energy on impact.

        **Velocity** — How fast it's moving. Energy grows with the **square** of speed!

        **Kinetic energy** — Energy of motion. We convert it to **megatons of TNT** to help compare sizes.

        **Crater** — A bowl-shaped hole. Our estimate is simplified for learning.

        **Deflection (Δv)** — A tiny push far in advance can move an asteroid enough to miss Earth.

        **Why simplified?** — Real scientists use more complex models (airbursts, fragmentation, terrain, oceans).
        Our goal is to **learn the ideas** first—then you can add advanced physics.
        """)
    st.markdown("### " + (t("app.where_to_extend") if t("app.where_to_extend", default=None) else t("app_extra.where_to_extend") or "Where to extend (hackathon tasks)"))
    st.markdown(t("app.where_to_extend_bullet_usgs") or "- **USGS overlays:** add coastal elevation/tsunami hazard layers via tiled map sources.")
    st.markdown(t("app.where_to_extend_bullet_population") or "- **Population exposure:** add a layer with night lights or population to illustrate risk.")
    st.markdown(t("app.where_to_extend_bullet_models") or "- **Better crater/blast models:** swap in published scaling relations and atmosphere effects.")
    st.markdown(t("app.where_to_extend_bullet_orbit") or "- **Orbit view:** a 3D Three.js canvas for the Sun–Earth–asteroid geometry (or use Plotly 3D).")

st.sidebar.title(t("sidebar.about_title") or "About this MVP")
if defaults_meta:
    default_name = (
        defaults_meta.get("sbdb_fullname")
        or defaults_meta.get("neo_name")
        or DEFAULT_SBDB_ID
    )
    provenance = defaults_meta.get("provenance", "NeoWs/SBDB")
    taxonomy = defaults_meta.get("taxonomy")
    density_val = defaults_meta.get("density")
    st.sidebar.markdown(
        "**" + t("sidebar.default_object", name=default_name) + "**<br>"
        + t("sidebar.source", src=provenance)
        + "<br>" + t("sidebar.spectral", tax=taxonomy or "n/a")
        + "<br>" + t("sidebar.density", dens=density_val or 0.0),
        unsafe_allow_html=True,
    )
st.sidebar.info(t("sidebar.disclaimer") or "This is an educational demo. Numbers are approximate. For real decisions, consult official models and data.")

control_state_snapshot = gather_ui_control_state()
print_startup_data_trace(defaults_meta or {}, control_state_snapshot)