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

from src.i18n import t, set_lang, get_lang, get_language_label, AVAILABLE_LANGS

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

# NOTE: These casualty rates are NOT from the PAIR white paper.
# Mathias et al. (2017) uses "affected population" (everyone within 4-psi damage radius)
# without modeling specific casualty rates. These rates are kept for legacy/educational purposes
# but should NOT be used in PAIR-compliant calculations.
SYNTHETIC_CASUALTY_RATE = {
    "severe": 0.35,
    "moderate": 0.1,
    "light": 0.02,
}

# Friendly city presets (lat, lon)
CITY_PRESETS = {
    "— choose a place —": (None, None),
    "Kaohsiung City, Taiwan": (22.6273, 120.3014),
    "Mulhouse, France": (47.7508, 7.3359),
    "Mexicali, Mexico": (32.6245, -115.4523),
    "Québec City, Canada": (46.8139, -71.2080),
    "Bangkok, Thailand": (13.7563, 100.5018),
}

# Approximate ring densities (people/km²) for educational use.
# Severe = inner blast zone, Moderate = suburban radius, Light = outer affected ring
CITY_RING_DENSITY = {
    "Kaohsiung City, Taiwan": {"severe": 9000, "moderate": 4000, "light": 1500},
    "Mulhouse, France":       {"severe": 4000, "moderate": 1500, "light": 600},
    "Mexicali, Mexico":       {"severe": 5000, "moderate": 2000, "light": 800},
    "Québec City, Canada":    {"severe": 4500, "moderate": 1800, "light": 700},
    "Bangkok, Thailand":      {"severe": 12000, "moderate": 5000, "light": 2000},
}

# Fallback for non-preset coordinates (rural/suburban)
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


def _console_print(*args: Any, **kwargs: Any) -> None:
    """Print to stdout but ignore BrokenPipe errors (Streamlit teardown)."""

    try:
        print(*args, **kwargs)
    except BrokenPipeError:  # pragma: no cover - depends on runtime teardown
        pass

def print_startup_data_trace(
    defaults: Dict[str, Any], 
    control_state: Optional[Dict[str, Any]] = None,
) -> None:
    global _DEFAULTS_TRACE_EMITTED
    if _DEFAULTS_TRACE_EMITTED:
        return
    _DEFAULTS_TRACE_EMITTED = True

    _console_print("\n=== Dashboard data provenance snapshot ===")
    neo_label = defaults.get("neo_name") or defaults.get("neo_designation") or DEFAULT_NEO_ID
    sbdb_label = defaults.get("sbdb_fullname") or DEFAULT_SBDB_ID
    provenance = defaults.get("provenance") or defaults.get("source", "synthetic fallback")
    if defaults.get("source_error"):
        _console_print(f"[warn] data fetch issues: {defaults['source_error']}")
    _console_print(f"NeoWs target : {neo_label} (ID {DEFAULT_NEO_ID})")
    _console_print(f"SBDB target  : {sbdb_label}")
    _console_print(f"Density source: {provenance}")

    field_sources = defaults.get("field_sources") or {}
    if field_sources:
        _console_print("Field-level sources:")
        for key in sorted(field_sources):
            _console_print(f"  - {key}: {field_sources[key]}")

    try:
        from tests.run_data_source_check import PARAMETERS as WHITEPAPER_PARAMS, get_nested
    except ImportError:
        _console_print("[info] Parameter comparison table unavailable (tests package not importable).")
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
        _console_print(f"[warn] NeoWs payload unavailable: {neo_error}")
    if sbdb_error:
        _console_print(f"[warn] SBDB payload unavailable: {sbdb_error}")

    selected_approach = defaults.get("close_approach_snapshot") or {}
    selected_velocity = defaults.get("velocity_km_s")
    selected_miss = defaults.get("miss_distance_km")

    header = f"{ 'Parameter':35} { 'NeoWs value':>18} { 'SBDB value':>18} { 'Dashboard':>18} { 'Source':>18}"
    _console_print("\n" + header)
    _console_print("-" * len(header))

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

        _console_print(
            f"{spec.label:35} "
            f"{_fmt(neo_value):>18} "
            f"{_fmt(sbdb_value):>18} "
            f"{_fmt(dashboard_value):>18} "
            f"{_fmt(source_label):>18}"
        )

    if control_state:
        _console_print("\nUI controls at startup:")
        for label, value in control_state.items():
            _console_print(f"  - {label}: {value}")


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
        altitude -= min(8.0, 6.0 * math.log10(max(diameter_m, 10.0) / 50.0))

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
    # Steeper logistic: midpoint ~15 km, scale ~2.5 km
    midpoint = 15.0
    scale = 2.5
    fraction = 1.0 / (1.0 + math.exp((burst_altitude_km - midpoint) / scale))
    # Guard against underflow/overflow and clamp
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

def final_crater_diameter_m(
    diameter_m: float,
    velocity_km_s: float,
    density_kg_m3: float,
    angle_deg: float,
    burst_altitude_km: float,
) -> float:
    """
    Compute a final crater diameter even if breakup occurs, by scaling with the
    ground-coupled energy fraction. Educational, not mission-grade.

    Rules of thumb:
    - If <5% of the kinetic energy reaches the ground → true airburst → no crater.
    - Otherwise, scale the transient crater by gf^(1/3) (crater size ~ E^(~1/3) scaling)
      and then apply a simple transient->final factor.
    """
    # Fraction of the entry kinetic energy that reaches the ground
    gf = ground_energy_fraction(burst_altitude_km)  # 0..1

    # True airburst: almost nothing reaches the ground
    if gf < AIRBURST_THRESHOLD:
        return 0.0

    # Compute a no-fragmentation transient crater
    D_t = transient_crater_diameter_m(diameter_m, velocity_km_s, density_kg_m3, angle_deg)
    if D_t <= 0:
        return 0.0

    # Reduce crater size by the fraction of energy that survives to the ground.
    # Crater diameter scales roughly with E^(1/3) in the gravity regime.
    D_t_eff = D_t * (gf ** (1.0 / 3.0))

    # If extremely small after scaling, treat as no excavated crater
    if D_t_eff < 50.0:  # meters; guard against tiny/ambiguous bowls
        return 0.0

    # Simple transient->final conversion
    return 1.25 * D_t_eff

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


def apply_deflection(lat, lon, delta_v_mm_s, lead_years, inbound_bearing_deg=0.0):
    """
    Educational deflection model with simple 'phase amplification':
    - A small Δv applied early shifts the encounter phase each orbit.
    - We approximate this leverage with a constant amplification factor (~3),
      representing how along-track timing errors accumulate by encounter time.
    - Lead time is in years.
    - If the amplified along-track displacement exceeds ~1 Earth radius, treat as a miss.

    Returns:
        (new_lat, new_lon, shift_km) if it still impacts,
        or (None, None, miss_distance_km) if it misses Earth.
    """
    dv = delta_v_mm_s / 1000.0          # mm/s -> m/s
    t = lead_years * 365.25 * 86400.0   # years -> s

    AMPLIFICATION = 3.0                 # tunable 2–5 for classroom intuition
    along_m = dv * t * AMPLIFICATION    # amplified along-track displacement (m)
    along_km = along_m / 1000.0

    miss_threshold_km = EARTH_RADIUS_KM * 1.05  # ~1 Earth radius + small margin

    if abs(along_km) >= miss_threshold_km:
        # Miss: return None coords and the (signed) miss distance in km
        return None, None, along_km

    # Still hits: shift ground impact point for visualization (cap so it stays readable)
    plotted_shift_km = max(0.0, min(abs(along_km), 5000.0))
    new_lat, new_lon = haversine_offset(lat, lon, plotted_shift_km, (inbound_bearing_deg + 180) % 360)
    return new_lat, new_lon, plotted_shift_km


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
    q = st.query_params.get("lang", None)
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

st.title("🛰️ Impactor-2025: Learn & Simulate")

# Initialize session defaults
ensure_session_defaults()
defaults_meta = st.session_state.get("defaults_metadata", {})

# Fetch today's NEOs
neo_df = fetch_today_neos()

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
        diameter_m = st.slider(
            "Asteroid diameter (m)",
            min_value=10,
            max_value=2000,
            value=diameter_default,
            step=10,
            key="diameter_m",
        )
    with primary_cols[1]:
        velocity_default = float(st.session_state.get("velocity_km_s", 18.0))
        velocity = st.slider(
            "Velocity at impact (km/s)",
            min_value=5.0,
            max_value=70.0,
            value=velocity_default,
            step=0.5,
            key="velocity_km_s",
        )
    with primary_cols[2]:
        angle_default = int(st.session_state.get("angle_deg", 45))
        angle = st.slider(
            "Impact angle (°)",
            min_value=10,
            max_value=90,
            value=angle_default,
            step=1,
            key="angle_deg",
        )
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
        density = st.slider(
            "Bulk density (kg/m³)",
            300,
            9000,
            value=int(preset_density),
            step=50,
            key=f"density_{material}",
        )
    with secondary_cols[1]:
        strength_mpa = st.slider(
            "Bulk compressive strength (MPa)",
            0.1,
            300.0,
            value=float(preset_strength),
            step=0.1,
            key=f"strength_{material}",
        )

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
    burst_alt_km = max(0.0, breakup_alt_km)
    ground_fraction = ground_energy_fraction(burst_alt_km)
    crater_m = final_crater_diameter_m(diameter_m, velocity, density, angle, burst_alt_km)
    # If very little energy reaches the ground, force airburst in the UI
    if ground_fraction < AIRBURST_THRESHOLD:
        crater_m = 0.0
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
    # PAIR uses 4-psi damage radius as the primary metric (white paper Section 2.3, line 236)
    # 4-psi is a full circle, not a ring, so we need to calculate the total area within r_mod
    pair_damage_area_km2 = math.pi * max(r_mod, 0.0) ** 2
    # Use the city/fallback density for the 4-psi region
    if preset in CITY_RING_DENSITY:
        pair_density = CITY_RING_DENSITY[preset].get("moderate", DEFAULT_RING_DENSITY["moderate"])
    else:
        pair_density = DEFAULT_RING_DENSITY["moderate"]
    pair_affected_population_4psi = pair_damage_area_km2 * pair_density

    cols3[0].metric("Population in 4-psi zone", f"{pair_affected_population_4psi:,.0f}")
    cols3[1].metric("Total exposed (all rings)", f"{exposure.get('total', 0.0):,.0f}")
    cols3[2].metric("Seismic Mw", f"{Mw:.1f}" if Mw is not None else "n/a")

    with st.expander("Damage assessment details"):
        # PAIR uses 4-psi as the primary damage threshold (full circle, not a ring)
        pair_damage_radius_km = r_mod  # 4-psi radius
        pair_damage_area_km2 = math.pi * max(pair_damage_radius_km, 0.0) ** 2
        # Population = full area within 4-psi × density
        pair_affected_pop = pair_damage_area_km2 * pair_density

        st.metric("Damage radius (4-psi)", f"{pair_damage_radius_km:.2f} km")
        st.metric("Damage area (full circle)", f"{pair_damage_area_km2:,.1f} km²")
        st.metric("Population", f"{pair_affected_pop:,.0f} people")

        # Show all three rings for reference/visualization
        st.markdown("**Blast overpressure rings:**")
        exposure_rows = []
        for ring, radius, psi in [("severe", r_severe, 12), ("moderate", r_mod, 4), ("light", r_light, 1)]:
            area = math.pi * max(radius, 0.0) ** 2
            pop = exposure.get(ring, 0.0)
            exposure_rows.append(
                {
                    "ring": ring,
                    "overpressure_psi": psi,
                    "radius_km": radius,
                    "area_km2": area,
                    "population": pop,
                }
            )
        st.dataframe(pd.DataFrame(exposure_rows))

    with st.expander("Probabilistic scenarios"):
        samples = st.slider("Number of Monte Carlo samples", 100, 1000, 300, step=50)
        if st.button("Run simulation", key="run_pair_button"):
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

    # Map visualization with crater and blast damage zones
    st.markdown("### Impact Visualization Map")

    # Map visualization with concentric circles
    st.session_state["impact_scenario"] = {
        "lat": lat, "lon": lon,
        "diameter_m": diameter_m, "density": density, "velocity": velocity,
        "angle": angle, "strength_mpa": strength_mpa,
        "E_mt": E_mt, "E_ground_mt": E_ground_mt,
        "burst_alt_km": burst_alt_km, "crater_km": crater_km,
        "r_severe": r_severe, "r_mod": r_mod, "r_light": r_light, "Mw": Mw,
        "exposure": exposure, "casualties": exposure.get("casualties", 0.0),
        "ring_densities": current_ring_densities,  # NEW
    }

    st.markdown("### Map")
    view_state = pdk.ViewState(latitude=lat, longitude=lon, zoom=6, bearing=0, pitch=30)

    def circle_layer(radius_km, color, name=""):
        return pdk.Layer(
            "ScatterplotLayer",
            data=pd.DataFrame({"lat": [lat], "lon": [lon], "name": [name]}),
            get_position="[lon, lat]",
            get_radius=radius_km * 1000,
            radius_min_pixels=1,
            radius_max_pixels=10000,
            get_fill_color=color,
            pickable=True,
            stroked=False,
            filled=True,
        )

    rings = [
        # Blast damage zones (from outer to inner)
        circle_layer(r_light, [255, 165, 0, 60], f"1-psi zone: {r_light:.1f} km"),
        circle_layer(r_mod, [255, 0, 0, 80], f"4-psi zone (PAIR): {r_mod:.1f} km"),
        circle_layer(r_severe, [139, 0, 0, 120], f"12-psi zone: {r_severe:.1f} km"),
        # Crater
        pdk.Layer(
            "ScatterplotLayer",
            data=pd.DataFrame({"lat": [lat], "lon": [lon], "name": [f"Crater: {crater_km:.2f} km diameter"]}),
            get_position="[lon, lat]",
            get_radius=crater_km * 500,  # crater radius in meters
            get_fill_color=[50, 50, 50, 220],
            get_line_color=[255, 255, 255, 255],
            line_width_min_pixels=2,
            pickable=True,
        ),
        # Impact point marker
        pdk.Layer(
            "ScatterplotLayer",
            data=pd.DataFrame({"lat": [lat], "lon": [lon], "name": ["Impact Point"]}),
            get_position="[lon, lat]",
            get_radius=1000,
            get_fill_color=[255, 255, 0, 255],
            get_line_color=[255, 0, 0, 255],
            line_width_min_pixels=2,
            pickable=True,
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

    # Map legend
    st.markdown("""
    **Map Legend:**
    - 🟡 **Yellow dot**: Impact point
    - ⚫ **Dark circle**: Crater ({crater_km:.2f} km diameter)
    - 🔴 **Dark red zone**: 12-psi overpressure (severe damage, {r_severe:.1f} km radius)
    - 🔴 **Red zone**: 4-psi overpressure (PAIR damage threshold, {r_mod:.1f} km radius)
    - 🟠 **Orange zone**: 1-psi overpressure (light damage, {r_light:.1f} km radius)
    """.format(crater_km=crater_km, r_severe=r_severe, r_mod=r_mod, r_light=r_light))

    # --- PATCH 2: make NEOs usable + educational summaries ---
    with st.expander("Use today's NASA NEOs as inputs"):
        df = neo_df.copy()
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
            neo_options = ["— use custom values —"] + df["name"].tolist()
            choice = st.selectbox(
                "Pick an asteroid to simulate (what-if it hit):",
                options=neo_options,
                key="neo_choice"
            )

            if choice != "— use custom values —":
                row = df.loc[df["name"] == choice].iloc[0]


                # Automatically update sliders when selection changes
                diameter_value = row["diameter_avg_m"]
                if diameter_value is None or pd.isna(diameter_value):
                    diameter_value = st.session_state.get("diameter_m", 150)
                diameter_value = int(np.clip(diameter_value, 10, 2000))

                velocity_value = row["velocity_km_s"]
                if velocity_value is None or pd.isna(velocity_value):
                    velocity_value = st.session_state.get("velocity_km_s", 18.0)
                velocity_value = float(np.clip(velocity_value, 5.0, 70.0))

                # Check if values need updating
                if (st.session_state.get("diameter_m") != diameter_value or
                    st.session_state.get("velocity_km_s") != velocity_value):
                    override_values = {
                        "diameter_m": diameter_value,
                        "velocity_km_s": velocity_value,
                        "angle_deg": 45,
                    }
                    st.session_state["widget_overrides"] = override_values
                    st.rerun()

                # Preview metric
                mass = asteroid_mass_kg(row["diameter_avg_m"], density)
                E_mt_preview = tnt_megatons(kinetic_energy_joules(mass, row["velocity_km_s"]))
                st.metric("Impact energy", f"{E_mt_preview:,.2f} Mt TNT")

with defend_tab:
    st.subheader("Try a deflection strategy ✨")
    st.write("See how a small push far in advance could change an asteroid’s impact point and effects on Earth.")

    # --- Load stored parameters from Explore tab ---
    if "impact_scenario" not in st.session_state:
        st.warning("⚠️ Please first create an asteroid scenario in the **Explore** tab.")
        st.stop()

    scenario = st.session_state["impact_scenario"]

    base_lat = scenario["lat"]
    base_lon = scenario["lon"]
    diameter_m = scenario["diameter_m"]
    density = scenario["density"]
    velocity = scenario["velocity"]
    angle = scenario["angle"]
    strength_mpa = scenario["strength_mpa"]
    E_mt_base = scenario["E_mt"]
    E_ground_mt_base = scenario["E_ground_mt"]
    burst_alt_km_base = scenario["burst_alt_km"]
    crater_km_base = scenario["crater_km"]
    r12_base = scenario["r_severe"]
    r4_base = scenario["r_mod"]
    r1_base = scenario["r_light"]
    Mw_base = scenario["Mw"]
    casualties = scenario["casualties"]
    ring_densities = scenario.get("ring_densities", DEFAULT_RING_DENSITY)

    # --- Deflection controls ---
    c1, c2, c3 = st.columns(3)
    with c1:
        delta_v_mm_s = st.slider("Δv (mm/s)", 0.0, 50.0, 1.0, step=0.5)
    with c2:
        lead_years = st.slider("Lead time (years)", 0.0, 50.0, 5.0, step=0.5)
    with c3:
        inbound_bearing = st.slider(t("deflect.bearing_label") or "Inbound bearing (°)", 0, 359, 90, key="deflect_bearing")

    # --- Compute new impact point ---
    new_lat, new_lon, shift_km = apply_deflection(base_lat, base_lon, delta_v_mm_s, lead_years, inbound_bearing)

    if new_lat is None:
        st.success(f"✅ Deflection successful! The asteroid misses Earth by ~{abs(shift_km):,.0f} km at encounter.")
        st.info("Try smaller Δv or shorter lead time to find the minimum nudge that still avoids impact.")
        st.stop()
    else:
        st.warning(f"⚠️ Still on impact course. Ground shift ≈ {shift_km:.1f} km.")

    # --- Recompute outcome at new site (same asteroid physics) ---
    burst_alt_km = estimate_burst_altitude_km(velocity, strength_mpa, diameter_m, angle)
    burst_alt_km = max(0.0, burst_alt_km)
    ground_fraction = ground_energy_fraction(burst_alt_km)
    
    m = asteroid_mass_kg(diameter_m, density)
    E_j = kinetic_energy_joules(m, velocity)
    E_mt = tnt_megatons(E_j)
    E_ground_mt = E_mt * ground_fraction
    crater_m = final_crater_diameter_m(diameter_m, velocity, density, angle, burst_alt_km)
    if ground_fraction < AIRBURST_THRESHOLD:
        crater_m = 0.0
    crater_km = crater_m / 1000.0

    r12 = blast_overpressure_radius_km(E_ground_mt, burst_alt_km, 12.0)
    r4 = blast_overpressure_radius_km(E_ground_mt, burst_alt_km, 4.0)
    r1 = blast_overpressure_radius_km(E_ground_mt, burst_alt_km, 1.0)
    Mw = seismic_moment_magnitude(E_j, ground_fraction=ground_fraction)

    exposure = estimate_population_impacts(r12, r4, r1, ring_densities=ring_densities)
    new_casualties = exposure.get("casualties", 0.0)
    casualty_diff = new_casualties - casualties

    # --- Comparison Table ---
    comp = pd.DataFrame({
        "Parameter": [
            "Impact latitude", "Impact longitude",
            "Breakup altitude (km)",
            "Kinetic energy (Mt TNT)", "Ground-coupled energy (Mt TNT)",
            "Crater diameter (km)",
            "Severe radius (12 psi, km)", "Moderate radius (4 psi, km)", "Light radius (1 psi, km)",
            "Seismic magnitude (Mw)",
            "Estimated casualties",
        ],
        "Before (original)": [
            f"{base_lat:.3f}", f"{base_lon:.3f}",
            f"{burst_alt_km_base:.1f}",
            f"{E_mt_base:.2f}", f"{E_ground_mt_base:.2f}",
            f"{crater_km_base:.2f}" if crater_km_base > 0 else "Airburst",
            f"{r12_base:.2f}", f"{r4_base:.2f}", f"{r1_base:.2f}",
            f"{Mw_base:.1f}" if Mw_base else "n/a",
            f"{casualties:.1f}" if casualties else "n/a",
        ],
        "After (deflected)": [
            f"{new_lat:.3f}", f"{new_lon:.3f}",
            f"{burst_alt_km:.1f}",
            f"{E_mt:.2f}", f"{E_ground_mt:.2f}",
            f"{crater_km:.2f}" if crater_km > 0 else "Airburst",
            f"{r12:.2f}", f"{r4:.2f}", f"{r1:.2f}",
            f"{Mw:.1f}" if Mw else "n/a",
            f"{new_casualties:.1f}" if new_casualties else "n/a",
        ]
    })
    st.markdown("### 🌍 Before vs After Deflection")
    st.dataframe(comp, use_container_width=True)

    # --- Map visualization (original + deflected impact) ---
    view_state = pdk.ViewState(latitude=base_lat, longitude=base_lon, zoom=5, pitch=30)

    def circle_layer(lat, lon, radius_km, color, name):
        return pdk.Layer(
            "ScatterplotLayer",
            data=pd.DataFrame({"lat": [lat], "lon": [lon]}),
            get_position="[lon, lat]",
            get_radius=radius_km * 1000,
            get_fill_color=color,
            stroked=False,
            pickable=False,
            filled=True,
            opacity=0.25,
            id=name,
        )

    layers = [
        # Original (red)
        circle_layer(base_lat, base_lon, r1_base, [255, 165, 0, 80], "light_before"),
        circle_layer(base_lat, base_lon, r4_base, [255, 0, 0, 100], "mod_before"),
        circle_layer(base_lat, base_lon, r12_base, [139, 0, 0, 150], "sev_before"),
        # Deflected (green)
        circle_layer(new_lat, new_lon, r1, [0, 255, 0, 60], "light_after"),
        circle_layer(new_lat, new_lon, r4, [0, 200, 100, 80], "mod_after"),
        circle_layer(new_lat, new_lon, r12, [0, 100, 0, 120], "sev_after"),
        # Line connecting them
        pdk.Layer(
            "LineLayer",
            data=pd.DataFrame({"source_lon": [base_lon], "source_lat": [base_lat],
                               "target_lon": [new_lon], "target_lat": [new_lat]}),
            get_source_position="[source_lon, source_lat]",
            get_target_position="[target_lon, target_lat]",
            get_width=3,
            get_color=[0, 200, 255, 200],
        ),
    ]

    if MAPBOX_TOKEN:
        deck = pdk.Deck(map_style="mapbox://styles/mapbox/dark-v11", initial_view_state=view_state, layers=layers)
    else:
        basemap = pdk.Layer(
            "TileLayer",
            data="https://tile.openstreetmap.org/{z}/{x}/{y}.png",
            min_zoom=0, max_zoom=19, tile_size=256,
        )
        deck = pdk.Deck(initial_view_state=view_state, layers=[basemap, *layers])

    st.pydeck_chart(deck)
    st.caption("Red = original impact | Green = deflected impact. Adjust Δv and lead time to see how even small pushes can change where — and how severely — an asteroid hits Earth.")

with learn_tab:
    import streamlit as st

    st.subheader("Asteroid Impact Glossary")




    # --------------------------
    # Your glossary JSON (pasted as-is)
    # --------------------------
    GLOSSARY_RAW = {
      "Airburst": {
        "term": "Airburst",
        "Definition": "An explosion in the atmosphere caused when an incoming body (asteroid or meteoroid) breaks apart under aerodynamic pressure before reaching the ground.",
        "Example": "The Chelyabinsk event in 2013 was an airburst: the object exploded at about 30 km altitude, causing a shockwave and damage on the ground.",
        "Fun fact": "An airburst can sometimes cause more ground-level damage than a small crater impact because the shockwave spreads over a wide area.",
        "URL": "https://assets.iflscience.com/assets/articleNo/72849/aImg/74031/touchdown-l.webp"
      },
      "Asteroid": {
        "term": "Asteroid",
        "Definition": "A rocky (or metallic) body orbiting the Sun, typically found in the asteroid belt between Mars and Jupiter; some cross Earth's orbit and become Near-Earth Objects (NEOs).",
        "Example": "The asteroid 99942 Apophis is a known Near-Earth Object under study for possible Earth close approaches.",
        "Fun fact": "The largest known asteroid, Ceres, is also classified as a dwarf planet and has a diameter of ~940 km.",
        "URL": "https://hips.hearstapps.com/hmg-prod/images/rock-on-starry-background-royalty-free-image-1645118962.jpg?resize=1200:*"
      },
      "Breakup altitude": {
        "term": "Breakup altitude",
        "Definition": "The height above ground where aerodynamic pressure (dynamic pressure) exceeds the structural strength of the body, causing it to fragment.",
        "Example": "If a stony asteroid breaks up at 20 km altitude, much of its energy may be dissipated before ground impact.",
        "Fun fact": "Breakup altitude depends not only on speed, but also on internal cracks and composition; a “rubble pile” asteroid breaks more easily.",
        "URL": "https://www.spacesafetymagazine.com/wp-content/uploads/2014/05/reentry-breakup.jpg"
      },
      "Crater": {
        "term": "Crater",
        "Definition": "The depression or cavity formed on the surface when an impactor (asteroid) strikes the ground, displacing material and ejecta.",
        "Example": "Meteor Crater in Arizona is about 1.2 km across and ~170 m deep, created about 50,000 years ago.",
        "Fun fact": "Some craters show ‘central peaks’—mountainous rebounds—if the impact was energetic enough to momentarily behave like a fluid.",
        "URL": "https://upload.wikimedia.org/wikipedia/commons/f/fd/Meteor_Crater_-_Arizona.jpg"
      },
      "Deflection (Δv)": {
        "term": "Deflection (Δv)",
        "Definition": "A small change in an asteroid’s velocity (Δv) applied early in its orbit to alter its path so it misses Earth.",
        "Example": "Adding a Δv of just a few mm/s years in advance can cause the asteroid to drift enough to avoid Earth.",
        "Fun fact": "In gravitational dynamics, a tiny push applied early often yields far greater effect than a large push applied late.",
        "URL": "https://qph.cf2.quoracdn.net/main-qimg-45551dbb9fbad51a769e897496294fd1"
      },
      "Diameter (D)": {
        "term": "Diameter (D)",
        "Definition": "The straight-line width of the asteroid (i.e. the maximum cross-sectional distance). Larger diameter generally means more mass and higher potential energy.",
        "Example": "An asteroid with diameter 100 m has vastly lower mass than one 1 km in diameter (volume scales with the cube).",
        "Fun fact": "If you lined up a 100 m asteroid next to a football stadium, it’d stretch across the entire field and beyond!",
        "URL": "https://upload.wikimedia.org/wikipedia/commons/7/7c/Circle_diameter.svg"
      },
      "Density (ρ)": {
        "term": "Density (ρ)",
        "Definition": "Mass per unit volume (ρ = m/V). For asteroids, typical values are ~3000 kg/m³ (stony) or ~7800 kg/m³ (iron).",
        "Example": "If an asteroid has volume 1 × 10⁶ m³ and density 3000 kg/m³, its mass is 3 × 10⁹ kg.",
        "Fun fact": "Some asteroids are “rubble piles”—collections of rocks loosely bound—so their bulk density is surprisingly low due to internal voids.",
        "URL": "https://blogger.googleusercontent.com/img/b/R29vZ2xl/AVvXsEgGiaxx7peo4zh1eh7uNmH5OkmdOMAXxsEu_kvjAyLy_ZGMr8sT4OPVBNnx0HRVCES9JG7liEhwemzc72lxxB5sKUTBFYly28RvBqWIEUfLAORYse2cRbLKPOOAlg1YrQ8txpXoqWHBBV-r/s1600/Chart+showing+types+of+asteroids+and+meteorites.jpg"
      },
      "Impact angle (θ)": {
        "term": "Impact angle (θ)",
        "Definition": "The angle between the asteroid’s trajectory and the local horizontal (ground). Shallow angles (e.g. < 15°) tend to produce long paths and airbursts; steep angles (near 90°) produce deep craters.",
        "Example": "An impact at θ = 30° will spread energy over a larger in-plane area than one at θ = 90°.",
        "Fun fact": "Most celestial impacts prefer angles around 45° because vertical or extreme shallow trajectories are statistically rare.",
        "URL": "https://i.sstatic.net/KRBHo.jpg"
      },
      "Kinetic energy (Eₖ)": {
        "term": "Kinetic energy (Eₖ)",
        "Definition": "The energy of the asteroid in motion, which is converted (upon impact or airburst) into heat, shock, light etc. Given by the equation Eₖ = ½ m v².",
        "Example": "If m = 1 × 10⁹ kg and v = 20,000 m/s, then Eₖ = ½ × 1×10⁹ × (2×10⁴)² = 2×10¹⁷ J (≈ 48 megatons of TNT).",
        "Fun fact": "Since Eₖ ∝ v², doubling speed quadruples energy — that’s why high-speed impacts are disproportionately destructive.",
        "URL": "https://img.jagranjosh.com/images/2024/July/1072024/kinetic-energy-definition-formula-derivation-types-examples-and-calculations.webp"
      },
      "Mass (m)": {
        "term": "Mass (m)",
        "Definition": "The amount of matter in the asteroid (m = density × volume). It determines momentum and, together with velocity, the kinetic energy.",
        "Example": "If ρ = 3000 kg/m³ and diameter = 100 m, volume ≈ (4/3)π(50 m)³ ≈ 5.24×10⁵ m³, so mass ≈ 1.57×10⁹ kg.",
        "Fun fact": "Mass is one of the hardest properties to measure remotely — scientists often infer it by observing how the asteroid perturbs nearby spacecraft or other bodies.",
        "URL": "https://energywavetheory.com/wp-content/uploads/2018/05/EarthMass.jpg"
      },
      "Monte Carlo simulation": {
        "term": "Monte Carlo simulation",
        "Definition": "A method that runs many randomized trials (varying size, velocity, angle, strength) to estimate a probabilistic risk distribution of impacts.",
        "Example": "We can simulate 10,000 possible asteroid paths and record how many hit Earth vs miss, creating a risk histogram.",
        "Fun fact": "The name ‘Monte Carlo’ comes from the casino in Monaco, because the method relies on randomness — like gambling draws.",
        "URL": "https://datascience.eu/wp-content/uploads/2020/03/monte_carlo_price_1-636x310-2.png"
      },
      "Overpressure radius": {
        "term": "Overpressure radius",
        "Definition": "The distance from the impact or explosion center within which the blast pressure exceeds certain thresholds (e.g. 12 psi, 4 psi, 1 psi) causing severe, moderate, or light damage.",
        "Example": "If overpressure of 4 psi is reached at 5 km, buildings within that radius may suffer significant damage.",
        "Fun fact": "Even a relatively small asteroid can produce overpressure over tens of kilometers — far larger than its actual size.",
        "URL": "https://static.wixstatic.com/media/603cad_d8d2a496375c41eca09a5368d6cf75b8~mv2.jpg/v1/fill/w_800,h_528,al_c,q_85,enc_avif,quality_auto/603cad_d8d2a496375c41eca09a5368d6cf75b8~mv2.jpg"
      },
      "Velocity (v)": {
        "term": "Velocity (v)",
        "Definition": "The speed of the asteroid relative to Earth (or the atmosphere) at impact. Since Eₖ ∝ v², velocity strongly impacts energy release.",
        "Example": "Typical asteroid entry speeds are 11–30 km/s; if v = 20,000 m/s and mass = 1×10⁹ kg, Eₖ = ~2×10¹⁷ J.",
        "Fun fact": "Even small increases in velocity yield huge gains in energy — a 10% speed boost gives ~21% more energy (since energy ∝ v²).",
        "URL": "https://media.hswstatic.com/eyJidWNrZXQiOiJjb250ZW50Lmhzd3N0YXRpYy5jb20iLCJrZXkiOiJnaWZcL3ZlbG9jaXR5LWZvcm11bGEtbmV3MS5qcGciLCJlZGl0cyI6eyJyZXNpemUiOnsid2lkdGgiOjgyOH0sInRvRm9ybWF0IjoiYXZpZiJ9fQ=="
      }
    }

    # --------------------------
    # Helpers
    # --------------------------
    def slugify(s: str) -> str:
        return "".join(ch.lower() if ch.isalnum() else "-" for ch in s).strip("-")

    def normalize_items(d: dict):
        """Convert your dict-of-dicts into a clean list of entries with consistent keys."""
        items = []
        for key, v in d.items():
            term = v.get("term", key)
            items.append({
                "term": term,
                "Definition": v.get("Definition", ""),
                "Example": v.get("Example", ""),
                "Fun fact": v.get("Fun fact", ""),
                "URL": v.get("URL", ""),
            })
        items.sort(key=lambda x: x["term"].lower())
        return items

    # Normalize and initialize state
    GLOSSARY = normalize_items(GLOSSARY_RAW)
    for item in GLOSSARY:
        st.session_state.setdefault(f"show-{slugify(item['term'])}", False)

    # --- Place the button here, after GLOSSARY is defined ---
    if st.button("📖 See all definitions"):
        # Check if all are open
        all_open = all(st.session_state[f"show-{slugify(item['term'])}"] for item in GLOSSARY)
        for item in GLOSSARY:
            st.session_state[f"show-{slugify(item['term'])}"] = not all_open
    # UI
    # --------------------------
    cols_per_row = 3
    gap = "large"

    for i, item in enumerate(GLOSSARY):
        if i % cols_per_row == 0:
            cols = st.columns(cols_per_row, gap=gap)

        with cols[i % cols_per_row]:
            with st.container(border=True):
                st.markdown(f"**{item['term']}**")
                key = f"show-{slugify(item['term'])}"
                with st.expander("Details", expanded=st.session_state[key]):
                    if item["Definition"]:
                        st.markdown(f"**Definition:** {item['Definition']}")
                    if item["Example"]:
                        st.markdown(f"**Example:** {item['Example']}")
                    if item["Fun fact"]:
                        st.markdown(f"**Fun fact:** {item['Fun fact']}")
                    url = item.get("URL") or ""
                    if isinstance(url, str) and url.strip():
                        st.image(url.strip(), use_container_width=True)

    st.markdown("## Important equations")

    IMPORTANT_EQUATIONS = {
      "Flight Dynamics Equations": {
        "term": "Flight Dynamics Equations",
        "Equation": [
          "dm/dt = -0.5 * ρ_air * v^3 * A * σ",
          "dv/dt = -0.5 * ρ_air * v^2 * A * C_D / m - g * sin(θ)",
          "dθ/dt = (v / (R_E + h) - g / v) * cos(θ)",
          "dh/dt = v * sin(θ)"
        ],
        "Output": "Describes how an asteroid's mass, velocity, flight path angle, and altitude change during its passage through Earth's atmosphere.",
        "Variable": [
          ["m", "Mass of the asteroid (kg)"],
          ["v", "Velocity of the asteroid (m/s)"],
          ["θ", "Flight path angle (radians)"],
          ["h", "Altitude (m)"],
          ["t", "Time (s)"],
          ["g", "Acceleration due to gravity (9.81 m/s²)"],
          ["ρ_air", "Local atmospheric density (kg/m³)"],
          ["R_E", "Radius of Earth (≈6371 km)"],
          ["A", "Cross-sectional area (m²)"],
          ["C_D", "Drag coefficient (dimensionless)"],
          ["σ", "Ablation coefficient (kg⁻¹ m²)"]
        ]
      },

      "Fragmentation and Dispersion": {
        "term": "Fragmentation and Dispersion",
        "Equation": [
          "ρ_air * v^2 > S",
          "S_child = S_parent * (m_parent / m_child)^a",
          "v_dispersion = v_cloud * sqrt((3.5 * ρ_air) / ρ_cloud)"
        ],
        "Output": "Determines when the asteroid breaks apart under aerodynamic pressure, how fragment strength scales, and how debris disperses through the atmosphere.",
        "Variable": [
          ["ρ_air", "Atmospheric density (kg/m³)"],
          ["v", "Asteroid velocity (m/s)"],
          ["S", "Aerodynamic strength (Pa)"],
          ["S_parent", "Strength of the parent fragment (Pa)"],
          ["S_child", "Strength of the child fragment (Pa)"],
          ["m_parent", "Mass of the parent fragment (kg)"],
          ["m_child", "Mass of the child fragment (kg)"],
          ["a", "Strength-scaling exponent (dimensionless)"],
          ["v_cloud", "Velocity of debris cloud (m/s)"],
          ["ρ_cloud", "Density of the cloud material (kg/m³)"]
        ]
      },

      "Blast Overpressure Damage": {
        "term": "Blast Overpressure Damage",
        "Equation": "R_ground = 2.09h - 0.449h²E^(-1/3) + 5.08E^(1/3)",
        "Output": "Estimates the ground radius affected by shockwave overpressure, indicating zones of structural damage or human injury.",
        "Variable": [
          ["R_ground", "Ground damage radius (km)"],
          ["h", "Burst altitude (km)"],
          ["E", "Impact energy (megaton TNT equivalent)"]
        ]
      },

      "Thermal Radiation Damage": {
        "term": "Thermal Radiation Damage",
        "Equation": [
          "r = sqrt((η * E) / (2π * Φ_j))",
          "R_ground = sqrt(r² - h²)"
        ],
        "Output": "Computes the distance where heat from an airburst causes third-degree burns and projects the affected radius onto the ground.",
        "Variable": [
          ["r", "Threshold radius from burst point (km)"],
          ["η", "Luminous efficiency (fraction of energy radiated as heat)"],
          ["E", "Impact energy (Joules or Mt)"],
          ["Φ_j", "Thermal exposure threshold (J/m²)"],
          ["R_ground", "Projected ground radius (km)"],
          ["h", "Burst altitude (km)"]
        ]
      },

      "Asteroid Diameter": {
        "term": "Asteroid Diameter",
        "Equation": "D = (1.326 × 10⁶) × 10^(-H/5) / sqrt(p_v)",
        "Output": "Estimates the asteroid’s physical diameter based on its observed brightness and surface reflectivity.",
        "Variable": [
          ["D", "Asteroid diameter (m)"],
          ["H", "Absolute magnitude (brightness)"],
          ["p_v", "Albedo (reflectivity coefficient)"]
        ]
      },

      "Entry Angle": {
        "term": "Entry Angle",
        "Equation": "θ° = (90° / π) * cos⁻¹(2U - 1)",
        "Output": "Determines the statistical entry angle of an asteroid as it enters the atmosphere, used in probabilistic simulations.",
        "Variable": [
          ["θ°", "Entry angle in degrees"],
          ["U", "Random number uniformly distributed between 0 and 1"]
        ]
      }
    }

    # Button to toggle all expanders
    def eq_slugify(s: str) -> str:
        return "".join(ch.lower() if ch.isalnum() else "-" for ch in s).strip("-")

    # Initialize session state for equations
    for eq in IMPORTANT_EQUATIONS.values():
        st.session_state.setdefault(f"show-eq-{eq_slugify(eq['term'])}", False)

    if st.button("📖 See all equations"):
        all_open = all(st.session_state[f"show-eq-{eq_slugify(eq['term'])}"] for eq in IMPORTANT_EQUATIONS.values())
        for eq in IMPORTANT_EQUATIONS.values():
            st.session_state[f"show-eq-{eq_slugify(eq['term'])}"] = not all_open

    # Render each equation in a container with expander
    for eq in IMPORTANT_EQUATIONS.values():
        key = f"show-eq-{eq_slugify(eq['term'])}"
        with st.container(border=True):
            st.markdown(f"**{eq['term']}**")
            with st.expander("Details", expanded=st.session_state[key]):
                # Equations
                if isinstance(eq["Equation"], list):
                    st.markdown("**Equation(s):**")
                    for e in eq["Equation"]:
                        st.latex(e)
                else:
                    st.markdown("**Equation:**")
                    st.latex(eq["Equation"])
                # Output
                st.markdown(f"**Output:** {eq['Output']}")
                # Variables
                st.markdown("**Variables:**")
                for var, desc in eq["Variable"]:
                    st.markdown(f"- `{var}`: {desc}")


    st.divider()

    st.markdown(
        """
        ### 🛠️ Future Classroom Add-Ons (Hackathon Ideas)
        - 🗺️ **USGS / NASA layers:** overlay coastlines, fault zones, or elevation for tsunami risk.
        - 👥 **Population exposure:** visualize how many people live near the impact area.
        - 🌡️ **Atmospheric effects:** add fireball brightness or shock-wave timing.
        - ☀️ **3D orbit view:** show the asteroid’s path around the Sun using Plotly 3D or Three.js.
        - 🎮 **Game mode:** give students missions — *“Save Earth with ≤ 1 mm/s Δv!”*
        """
    )

    st.caption("Educational mode: simplified for learning. Data and models inspired by NASA, ESA, and academic impact simulations.")



st.sidebar.title(t("sidebar.about_title") or "About this MVP")

# Language selector in sidebar
# Use the available languages from the i18n module
if AVAILABLE_LANGS:
    available_codes = AVAILABLE_LANGS
else:
    available_codes = ["en", "es", "fr", "zh-Hant"]

session_lang = st.session_state.get("lang")
current_lang = session_lang or get_lang()

if current_lang not in available_codes:
    current_lang = available_codes[0]
    st.session_state["lang"] = current_lang

# Ensure the active translator matches the session state
if current_lang != get_lang():
    try:
        set_lang(current_lang)
    except ValueError:
        current_lang = available_codes[0]
        st.session_state["lang"] = current_lang
        set_lang(current_lang)

# Find current index for the select box
try:
    current_idx = available_codes.index(current_lang)
except ValueError:
    current_idx = 0

selected_lang = st.sidebar.selectbox(
    "Language / Idioma / Langue / 語言",
    options=available_codes,
    format_func=lambda code: get_language_label(code, fallback=code),
    index=current_idx,
    key="language_selector"
)

# Update language if changed
if selected_lang != current_lang:
    set_lang(selected_lang)
    st.session_state["lang"] = selected_lang
    st.rerun()

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
