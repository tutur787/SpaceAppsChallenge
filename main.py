# app.py
import math
import requests
import numpy as np
import pandas as pd
import streamlit as st
from folium import Map, Circle, Marker, LayerControl, LatLngPopup
from streamlit_folium import st_folium

# ---------- Page config ----------
st.set_page_config(page_title="Impactor-2025 Simulator", layout="wide")
st.title("Impactor-2025: Interactive Impact & Mitigation Simulator (Python-only)")
st.caption("Educational, first-order estimates. Not for real-world hazard planning.")

# ---------- Sidebar inputs ----------
st.sidebar.header("Asteroid Parameters")
use_nasa = st.sidebar.checkbox("Fetch from NASA NEO API", value=False)
diam_m = st.sidebar.number_input("Diameter (m)", min_value=50.0, max_value=5000.0, value=300.0, step=10.0)
rho = st.sidebar.number_input("Density (kg/m³)", min_value=500.0, max_value=8000.0, value=3000.0, step=100.0)
v_kms = st.sidebar.number_input("Impact velocity (km/s)", min_value=5.0, max_value=72.0, value=20.0, step=0.5)
impact_angle = st.sidebar.slider("Impact angle from horizontal (°)", 5, 90, 45)
is_ocean = st.sidebar.checkbox("Ocean impact", value=False)

st.sidebar.header("Mitigation (Δv deflection)")
dv_mm_s = st.sidebar.slider("Δv (mm/s)", 0, 500, 0)
lead_time_days = st.sidebar.slider("Lead time before encounter (days)", 0, 3650, 0)

# ---------- NASA NEO (optional) ----------
def fetch_neo_by_date(date_str: str, api_key: str = "QUInDCVww4vIfJPWDayrXmbUN76wv9jTKGVZzola") -> pd.DataFrame:
    """Fetch near-Earth objects for a date and return a small DataFrame."""
    url = "https://api.nasa.gov/neo/rest/v1/feed"
    try:
        resp = requests.get(
            url,
            params={"start_date": date_str, "end_date": date_str, "api_key": api_key},
            timeout=20,
        )
        resp.raise_for_status()
        data = resp.json()
        neos = []
        for _, objs in data.get("near_earth_objects", {}).items():
            for o in objs:
                try:
                    est = o["estimated_diameter"]["meters"]
                    mean_d = 0.5 * (est["estimated_diameter_min"] + est["estimated_diameter_max"])
                    relv = float(o["close_approach_data"][0]["relative_velocity"]["kilometers_per_second"])
                    neos.append(
                        {
                            "name": o.get("name", "unknown"),
                            "diam_m": mean_d,
                            "v_kms": relv,
                        }
                    )
                except Exception:
                    continue
        return pd.DataFrame(neos)
    except Exception as e:
        st.warning(f"NASA API request failed: {e}")
        return pd.DataFrame()

if use_nasa:
    st.sidebar.write("Example NASA fetch (by date range):")
    date = st.sidebar.text_input("Date (YYYY-MM-DD)", "2025-10-04")
    api_key = st.sidebar.text_input("NASA API key", "DEMO_KEY")
    if st.sidebar.button("Fetch NEOs"):
        df = fetch_neo_by_date(date, api_key)
        if df.empty:
            st.warning("No NEOs returned. Using manual inputs.")
        else:
            st.dataframe(df, use_container_width=True)
            pick = st.selectbox("Apply NEO row:", df.index, key="neo_row_pick")
            if pick is not None and not df.empty:
                diam_m = float(df.loc[pick, "diam_m"])
                v_kms = float(df.loc[pick, "v_kms"])
                st.success(f"Applied NEO: diameter={diam_m:.1f} m, v={v_kms:.2f} km/s")

# ---------- Physics (first-order) ----------
R = diam_m / 2.0
volume = (4.0 / 3.0) * math.pi * R**3
mass = rho * volume  # kg
v = v_kms * 1000.0   # m/s
E = 0.5 * mass * v**2  # Joules
MT_TNT = E / 4.184e15  # megatons TNT (1 Mt = 4.184e15 J)

st.subheader("Impact energetics")
col1, col2, col3 = st.columns(3)
col1.metric("Mass (kg)", f"{mass:,.0f}")
col2.metric("Impact energy (J)", f"{E:,.3e}")
col3.metric("Impact energy (Mt TNT)", f"{MT_TNT:,.1f}")
st.caption("Energetics and damage zones are illustrative.")

# Very simple damage radii (illustrative placeholders)
# Scale rings loosely vs. energy^(1/3); coefficients are tuned for demos
severe_km = 2.0 * (MT_TNT ** (1 / 3.0)) / 10.0
moderate_km = 5.0 * (MT_TNT ** (1 / 3.0)) / 10.0
light_km = 10.0 * (MT_TNT ** (1 / 3.0)) / 10.0

# ---------- Map & impact selection ----------
st.subheader("Choose Impact Location")

# First map: world view for selecting a point
m_world = Map(location=[20, 0], zoom_start=2, tiles="OpenStreetMap")
m_world.add_child(LatLngPopup())  # visual feedback on click

# IMPORTANT: do not suppress return; give the widget a stable key
map_state = st_folium(m_world, height=450, key="impact_map_world")

# Get clicked lat/lon from widget state or session
latlon = None
if map_state and map_state.get("last_clicked"):
    latlon = (map_state["last_clicked"]["lat"], map_state["last_clicked"]["lng"])
elif st.session_state.get("latlon"):
    latlon = st.session_state["latlon"]

if latlon:
    st.session_state["latlon"] = latlon

    # Zoomed-in map centered on chosen impact point
    m_zoom = Map(location=latlon, zoom_start=6, tiles="OpenStreetMap")
    Marker(latlon, tooltip="Impact point").add_to(m_zoom)

    # Draw damage rings
    for r_km, color, label in [
        (severe_km, "red", "Severe"),
        (moderate_km, "orange", "Moderate"),
        (light_km, "blue", "Light"),
    ]:
        Circle(
            latlon,
            radius=max(r_km * 1000.0, 1.0),  # meters; guard tiny radii
            color=color,
            fill=True,
            fill_opacity=0.15,
            popup=f"{label} ~{r_km:.1f} km",
        ).add_to(m_zoom)

    # Simple ocean toggle cue (illustrative only)
    if is_ocean:
        Circle(
            latlon,
            radius=max(light_km * 1500.0, 1.0),
            color="green",
            fill=False,
            popup="Illustrative tsunami band",
        ).add_to(m_zoom)

    LayerControl().add_to(m_zoom)
    st_folium(m_zoom, height=500, key="impact_map_zoomed")
else:
    st.info("Click anywhere on the world map to set the impact point.")

# ---------- Deflection toy model ----------
st.subheader("Deflection outcome (toy)")
if dv_mm_s > 0 and lead_time_days > 0:
    dv = dv_mm_s / 1000.0  # m/s
    t = lead_time_days * 86400
    # Very simplified along-track displacement
    miss_km = (dv * t) / 1000.0
    st.write(f"Estimated along-track displacement at encounter: **~{miss_km:,.0f} km** (toy model).")
    if miss_km > 10000:
        st.success("Qualitatively: likely Earth miss in this toy scenario.")
    else:
        st.info("Qualitatively: may still be an impact without additional Δv or lead time.")
else:
    st.caption("Set Δv and lead time to explore a qualitative deflection outcome.")
