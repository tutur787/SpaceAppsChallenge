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


def estimate_burst_altitude_km(velocity_km_s: float, strength_mpa: float) -> float:
    """Invert the PAIR breakup criterion (Mathias et al. 2017 Eq. 2) with an exponential atmosphere."""
    if velocity_km_s <= 0 or strength_mpa <= 0:
        return 0.0
    velocity_m_s = velocity_km_s * 1000.0
    strength_pa = strength_mpa * 1e6
    rho_break = strength_pa / max(velocity_m_s ** 2, 1.0)
    if rho_break >= SEA_LEVEL_DENSITY:
        return 0.0  # survives to ground before meeting breakup criterion
    rho_break = max(rho_break, 1e-9)
    altitude = -SCALE_HEIGHT_KM * math.log(rho_break / SEA_LEVEL_DENSITY)
    return float(max(0.0, min(altitude, 80.0)))


def ground_energy_fraction(burst_altitude_km: float) -> float:
    """Approximate fraction of kinetic energy that couples to the ground (PAIR-inspired)."""
    if burst_altitude_km <= 0:
        return 1.0
    if burst_altitude_km >= 25.0:
        return 0.0
    return max(0.0, 1.0 - burst_altitude_km / 25.0)


def transient_crater_diameter_m(
    diameter_m: float,
    velocity_km_s: float,
    density_kg_m3: float,
    angle_deg: float,
    target_density_kg_m3: float = RHO_TARGET,
    gravity: float = G,
) -> float:
    if diameter_m <= 0 or velocity_km_s <= 0:
        return 0.0
    v = velocity_km_s * 1000.0
    angle_rad = math.radians(max(5.0, min(angle_deg, 90.0)))
    angle_factor = math.sin(angle_rad) ** (1 / 3)
    gravity_term = (gravity * diameter_m / (v ** 2)) ** -0.217
    density_term = (density_kg_m3 / target_density_kg_m3) ** 0.333
    return 1.161 * gravity_term * density_term * (diameter_m ** 0.783) * angle_factor


def final_crater_diameter_m(
    diameter_m: float,
    velocity_km_s: float,
    density_kg_m3: float,
    angle_deg: float,
    burst_altitude_km: float,
) -> float:
    energy_fraction = ground_energy_fraction(burst_altitude_km)
    if energy_fraction <= 0:
        return 0.0
    transient = transient_crater_diameter_m(diameter_m, velocity_km_s, density_kg_m3, angle_deg)
    if transient <= 0:
        return 0.0
    transient *= energy_fraction ** (1 / 3)
    if transient < 3200:
        final_crater = 1.25 * transient
    else:
        final_crater = 1.17 * (transient ** 1.13) / (G ** 0.13)
    return float(max(0.0, final_crater))


def blast_overpressure_radius_km(E_mt: float, burst_altitude_km: float, overpressure_psi: float = 4.0) -> float:
    if E_mt <= 0:
        return 0.0
    E13 = max(E_mt, 1e-6) ** (1 / 3)
    h = max(0.0, burst_altitude_km)
    base_radius = 2.09 * h - 0.449 * h * h / E13 + 5.08 * E13
    base_radius = max(0.0, base_radius)
    if overpressure_psi == 4.0:
        return base_radius
    scale = (4.0 / max(overpressure_psi, 1e-3)) ** (1 / 3)
    return base_radius * scale


def seismic_moment_magnitude(E_joules: float, coupling: float = DEFAULT_SEISMIC_COUPLING) -> Optional[float]:
    if E_joules <= 0 or coupling <= 0:
        return None
    seismic_energy = E_joules * coupling
    if seismic_energy <= 0:
        return None
    return (2.0 / 3.0) * (math.log10(seismic_energy) - 4.8)


def estimate_population_impacts(r_severe: float, r_moderate: float, r_light: float) -> dict[str, float]:
    radii = sorted([r_severe, r_moderate, r_light])
    if radii[-1] <= 0:
        return {"severe": 0.0, "moderate": 0.0, "light": 0.0}
    severe = max(r_severe, 0.0)
    moderate = max(r_moderate, severe)
    light = max(r_light, moderate)

    areas = {
        "severe": math.pi * severe ** 2,
        "moderate": math.pi * moderate ** 2 - math.pi * severe ** 2,
        "light": math.pi * light ** 2 - math.pi * moderate ** 2,
    }

    exposure = {}
    for ring, area in areas.items():
        density = SYNTHETIC_POP_DENSITY.get(ring, 0)
        exposure[ring] = max(0.0, area) * density
    exposure["total"] = sum(exposure.values())

    casualties = 0.0
    for ring, pop in exposure.items():
        if ring == "total":
            continue
        rate = SYNTHETIC_CASUALTY_RATE.get(ring, 0.0)
        casualties += pop * rate
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
        burst = estimate_burst_altitude_km(v, s)
        ground_frac = ground_energy_fraction(burst)
        crater = final_crater_diameter_m(d, v, rho, theta, burst) / 1000.0
        r4 = blast_overpressure_radius_km(E_mt, burst, overpressure_psi=4.0)
        r12 = blast_overpressure_radius_km(E_mt, burst, overpressure_psi=12.0)
        r1 = blast_overpressure_radius_km(E_mt, burst, overpressure_psi=1.0)
        exposure = estimate_population_impacts(r12, r4, r1)
        Mw = seismic_moment_magnitude(E_j)
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
st.caption("An educational dashboard to explore asteroid impacts, built for a hackathon.")

# Tabs
exp_tab, defend_tab, learn_tab = st.tabs(["Explore", "Defend Earth", "Learn"])

with exp_tab:
    st.subheader("Choose impact parameters")
    primary_cols = st.columns(4)
    with primary_cols[0]:
        diameter_m = st.slider("Asteroid diameter (m)", 10, 2000, 150, step=10)
    with primary_cols[1]:
        velocity = st.slider("Velocity at impact (km/s)", 5.0, 70.0, 18.0, step=0.5)
    with primary_cols[2]:
        angle = st.slider("Impact angle (°)", 10, 90, 45, step=1)
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

    # Computations
    m = asteroid_mass_kg(diameter_m, density)
    E_j = kinetic_energy_joules(m, velocity)
    E_mt = tnt_megatons(E_j)
    burst_alt_km = estimate_burst_altitude_km(velocity, strength_mpa)
    ground_fraction = ground_energy_fraction(burst_alt_km)
    crater_m = final_crater_diameter_m(diameter_m, velocity, density, angle, burst_alt_km)
    crater_km = crater_m / 1000.0
    r_mod = blast_overpressure_radius_km(E_mt, burst_alt_km, overpressure_psi=4.0)
    r_severe = blast_overpressure_radius_km(E_mt, burst_alt_km, overpressure_psi=12.0)
    r_light = blast_overpressure_radius_km(E_mt, burst_alt_km, overpressure_psi=1.0)
    exposure = estimate_population_impacts(r_severe, r_mod, r_light)
    Mw = seismic_moment_magnitude(E_j)

    st.markdown("### Results")
    cols = st.columns(4)
    cols[0].metric("Mass", f"{m:,.0f} kg")
    cols[1].metric("Kinetic energy", f"{E_mt:,.2f} Mt TNT")
    cols[2].metric("Breakup altitude", f"{burst_alt_km:.1f} km")
    cols[3].metric("Ground-coupled energy", f"{E_mt * ground_fraction:.2f} Mt")

    crater_display = f"{crater_km:.2f} km" if crater_km > 0.0 else "Airburst"
    cols2 = st.columns(4)
    cols2[0].metric("Crater diameter (final)", crater_display)
    cols2[1].metric("Severe radius (12 psi)", f"{r_severe:.2f} km")
    cols2[2].metric("Moderate radius (4 psi)", f"{r_mod:.2f} km")
    cols2[3].metric("Light radius (1 psi)", f"{r_light:.2f} km")

    cols3 = st.columns(3)
    cols3[0].metric("Population exposed", f"{exposure['total']:,.0f}")
    cols3[1].metric("Estimated casualties", f"{exposure['casualties']:,.0f}")
    cols3[2].metric("Seismic Mw", f"{Mw:.1f}" if Mw is not None else "n/a")

    st.caption("Population estimates rely on documented synthetic density rings; replace with real datasets when available.")

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

    with st.expander("Optional: fetch today's NEOs from NASA (for context)"):
        df = fetch_today_neos()
        if df.empty:
            st.info("Could not fetch NEOs right now (API rate limit or offline). You can still use the simulator.")
        else:
            st.dataframe(df)

with defend_tab:
    st.subheader("Try a deflection strategy ✨")
    st.write("Toy model: apply a small velocity change (Δv) some days before arrival and see how the nominal impact point shifts.")

    c1, c2, c3 = st.columns(3)
    with c1:
        delta_v_mm_s = st.slider("Δv (mm/s)", 0.0, 5.0, 0.5, step=0.1)
    with c2:
        lead_days = st.slider("Lead time (days)", 0, 3650, 365, step=30)
    with c3:
        inbound_bearing = st.slider("Inbound bearing (°)", 0, 359, 90)

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

st.sidebar.title("About this MVP")
st.sidebar.info(
    "This is an educational demo. Numbers are approximate. For real decisions, consult official models and data."
)
