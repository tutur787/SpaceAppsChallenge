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
from typing import Optional

import numpy as np
import pandas as pd
import requests
import streamlit as st
import pydeck as pdk

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

# Airburst decision: if less than this fraction of KE reaches the ground, treat as no crater
AIRBURST_THRESHOLD = 0.05

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
    diameter_m: float | None = None,
    angle_deg: float | None = None,   # NEW
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

def fetch_today_neos():
    key = 'QUInDCVww4vIfJPWDayrXmbUN76wv9jTKGVZzola'
    today = date.today().isoformat()
    url = f"https://api.nasa.gov/neo/rest/v1/feed?start_date={today}&end_date={today}&api_key={key}"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.json()
        neos = []
        for neo in data.get("near_earth_objects", {}).get(today, []):
            neos.append({
                "name": neo.get("name"),
                "est_diameter_min_m": neo["estimated_diameter"]["meters"]["estimated_diameter_min"],
                "est_diameter_max_m": neo["estimated_diameter"]["meters"]["estimated_diameter_max"],
                "hazardous": neo.get("is_potentially_hazardous_asteroid", False),
                "velocity_km_s": float(neo["close_approach_data"][0]["relative_velocity"]["kilometers_per_second"]) if neo.get("close_approach_data") else None,
                "miss_distance_km": float(neo["close_approach_data"][0]["miss_distance"]["kilometers"]) if neo.get("close_approach_data") else None,
            })
        return pd.DataFrame(neos)
    except Exception:
        return pd.DataFrame([])


# ----------------------------
# UI
# ----------------------------
st.set_page_config(page_title="Impactor-2025: Learn & Simulate", layout="wide")
st.title("🛰️ Impactor-2025: Learn & Simulate")
st.caption("An educational dashboard to explore asteroid impacts.")
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
        diameter_m = st.slider("Asteroid diameter (m)", **diameter_kwargs)
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
        velocity = st.slider("Velocity at impact (km/s)", **velocity_kwargs)
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
        angle = st.slider("Impact angle (°)", **angle_kwargs)
    with primary_cols[3]:
        material = st.selectbox("Material preset", list(MATERIAL_PRESETS.keys()), index=1)

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
        st.caption(f"Preset density for {material}: {preset_density} kg/m³")
    with secondary_cols[1]:
        strength_mpa = st.slider(
            "Bulk compressive strength (MPa)",
            0.1,
            300.0,
            value=float(preset_strength),
            step=0.1,
            key=f"strength_{material}",
        )
        st.caption("Adjust to emulate cohesive strength used in PAIR entry modeling.")

    st.markdown("**Where does it hit?** Pick a city or enter coordinates.")
    c1, c2, c3 = st.columns([2,1,1])
    with c1:
        preset = st.selectbox("City preset", list(CITY_PRESETS.keys()))
    if preset and CITY_PRESETS[preset][0] is not None:
        lat, lon = CITY_PRESETS[preset]
    else:
        with c2:
            lat = st.number_input("Latitude", value=29.7604, format="%.4f")
        with c3:
            lon = st.number_input("Longitude", value=-95.3698, format="%.4f")
    # Choose densities based on city selection
    if preset in CITY_RING_DENSITY:
        current_ring_densities = CITY_RING_DENSITY[preset]
    else:
        current_ring_densities = DEFAULT_RING_DENSITY

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
        samples = st.slider("Number of Monte Carlo samples", 100, 1000, 300, step=50)
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
    with st.expander("Use today's NASA NEOs as inputs"):
        df = fetch_today_neos()
        if df.empty:
            st.info("Could not fetch NEOs right now (API rate limit or offline). You can still use the simulator.")
        else:
            # Derive average diameter and lunar-distance multiples for context
            df = df.copy()
            df["diameter_avg_m"] = 0.5 * (df["est_diameter_min_m"] + df["est_diameter_max_m"])
            MOON_KM = 384_400.0
            df["miss_distance_moon_x"] = df["miss_distance_km"] / MOON_KM

            # Short, readable view
            view_cols = ["name", "diameter_avg_m", "velocity_km_s", "miss_distance_km", "miss_distance_moon_x", "hazardous"]
            st.dataframe(df[view_cols].rename(columns={
                "diameter_avg_m": "diameter_avg_m (m)",
                "miss_distance_km": "miss_distance (km)",
                "miss_distance_moon_x": "miss distance (× Moon)"
            }))

            # Pick one and push into the simulator
            choice = st.selectbox(
                "Pick an asteroid to simulate (what-if it hit):",
                options=df["name"].tolist()
            )
            row = df.loc[df["name"] == choice].iloc[0]

            st.caption(
                f"Selected **{row['name']}** — avg diameter ≈ {row['diameter_avg_m']:.1f} m, "
                f"speed ≈ {row['velocity_km_s']:.2f} km/s, "
                f"miss distance ≈ {row['miss_distance_km']:.0f} km (~{row['miss_distance_moon_x']:.1f}× Moon)."
            )

            colA, colB = st.columns(2)
            with colA:
                if st.button("Use this NEO in the simulator"):
                    diameter_value = row["diameter_avg_m"]
                    if diameter_value is None or pd.isna(diameter_value):
                        diameter_value = st.session_state["diameter_m"]
                    diameter_value = int(np.clip(diameter_value, 10, 2000))

                    velocity_value = row["velocity_km_s"]
                    if velocity_value is None or pd.isna(velocity_value):
                        velocity_value = st.session_state["velocity_km_s"]
                    velocity_value = float(np.clip(velocity_value, 5.0, 70.0))

                    override_values = {
                        "diameter_m": diameter_value,
                        "velocity_km_s": velocity_value,
                        # Reset to a representative entry angle for clarity when swapping asteroids
                        "angle_deg": 45,
                    }
                    st.session_state["widget_overrides"] = override_values
                    st.rerun()
            with colB:
                # A quick educational “scale” card
                mass = asteroid_mass_kg(row["diameter_avg_m"], density)
                E_mt_preview = tnt_megatons(kinetic_energy_joules(mass, row["velocity_km_s"]))
                st.metric("What-if energy (preview)", f"{E_mt_preview:,.2f} Mt TNT")
                st.caption("Preview assumes current density/material selection.")

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
        inbound_bearing = st.slider("Inbound bearing (°)", 0, 359, 90)

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
    st.header("🪐 Learn About Asteroids & Planetary Defense")

    st.markdown(
        """
        ### ☄️ What Are Asteroids?
        Asteroids are rocky or metallic leftovers from the early Solar System — like **building blocks** that never formed a planet.
        Most orbit safely between **Mars and Jupiter**, but some drift close to Earth.  
        These are called **Near-Earth Objects (NEOs)**.

        **Fun fact:** Scientists track over *35,000 NEOs*! NASA’s Planetary Defense Office keeps watch for us every day.
        """
    )

    st.divider()

    st.markdown(
        """
        ### ⚙️ Key Terms

        | Concept | What It Means | Why It Matters |
        |:---------|:---------------|:----------------|
        | **Diameter** | How wide the asteroid is (meters or kilometers). | Bigger asteroids carry **more mass** and release **more energy**. |
        | **Velocity** | How fast it's moving (km/s). | Energy increases with the **square of speed** — doubling speed means **4×** the energy! |
        | **Density (ρ)** | How heavy it is for its size. | Iron asteroids hit harder than icy ones. |
        | **Kinetic Energy** | Energy of motion: ½ × m × v² | Used to estimate damage, expressed in **megatons of TNT**. |
        | **Breakup Altitude** | Where it explodes in the air. | High breakup → airburst (like Chelyabinsk). Low breakup → crater (like Barringer). |
        | **Crater** | The bowl left by an impact. | Its size tells us about the energy released. |
        | **Deflection (Δv)** | A gentle nudge applied long before impact. | Even a few **millimeters per second** can make Earth miss the asteroid! |
        """
    )

    st.info(
        "💡 **Analogy:** Think of deflection like steering a car early — "
        "a tiny turn far away can mean a huge difference down the road!"
    )

    st.divider()

    st.markdown(
        """
        ### 🎓 For Students & Educators
        - **Elementary (Ages 8–11):**  
          Explore how size and speed change impact power. Ask: “What if it were as big as our school?”
        - **Middle School (Ages 12–14):**  
          Try graphing how **energy** changes with **velocity**. Notice it’s not linear!
        - **High School (Ages 15–18):**  
          Use the equations to estimate real events — *Chelyabinsk*, *Tunguska*, *Barringer*.  
          Discuss mitigation: **kinetic impactors**, **gravity tractors**, **laser ablation**, or **nuclear options**.
        """
    )

    st.markdown(
        """
        **Try it! 🔬**  
        - Create your own asteroid in the **Explore** tab.  
        - Use **Defend** to see how different Δv values change where it lands.  
        - Challenge: find the smallest Δv that moves it safely off Earth!
        """
    )

    st.divider()

    st.markdown(
        """
        ### 🧠 Why Our Model Is Simplified
        Real scientists include:
        - **Fragmentation physics** — asteroids can break apart like gravel.
        - **Atmospheric drag and ablation** — they lose mass and slow down.
        - **Topography & oceans** — mountains and water change crater and tsunami effects.
        - **Detailed scaling laws** — from NASA, USGS, and the *Earth Impact Effects Program*.

        Our dashboard focuses on **concepts**, not perfect precision — just like a flight simulator teaches before you fly a real plane. ✈️
        """
    )

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


st.sidebar.title("About this MVP")
st.sidebar.info(
    "This is an educational demo. Numbers are approximate. For real decisions, consult official models and data."
)
