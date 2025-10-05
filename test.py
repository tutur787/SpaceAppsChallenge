import math
import pandas as pd

# --- Constants ---
G = 9.81                    # m/s²
SEA_LEVEL_DENSITY = 1.225   # kg/m³
SCALE_HEIGHT = 8000.0       # m
TNT_J = 4.184e9             # J/ton TNT

# --- Helper functions ---
def kinetic_energy_mt(diameter_m, rho_i, v_kms):
    r = diameter_m / 2
    m = (4/3) * math.pi * r**3 * rho_i
    v = v_kms * 1000
    E_J = 0.5 * m * v**2
    return E_J / (TNT_J * 1e6)  # megatons TNT

def breakup_altitude_km(strength_mpa, rho_i, v_kms):
    """Approximate altitude where dynamic pressure = material strength"""
    strength_Pa = strength_mpa * 1e6
    v = v_kms * 1000
    rho_break = 2 * strength_Pa / (v**2)
    if rho_break > SEA_LEVEL_DENSITY:
        return 0.0  # hits ground before breaking
    h_m = -SCALE_HEIGHT * math.log(rho_break / SEA_LEVEL_DENSITY)
    return h_m / 1000

def crater_diameter_km(E_mt):
    """Rough scaling for simple craters"""
    return 0.01 * (E_mt)**0.3   # empirical (fits Barringer, Sikhote-Alin, etc.)

def seismic_magnitude(E_mt):
    """Approximate seismic equivalent magnitude"""
    E_J = E_mt * TNT_J * 1e6
    return 0.67 * math.log10(E_J) - 5.87

def classify_event(breakup_km):
    if breakup_km > 5:
        return "Airburst"
    elif breakup_km > 0.5:
        return "Low-altitude airburst"
    else:
        return "Surface impact"

# --- Run scenarios ---
scenarios = [
    ("Chelyabinsk (2013)", 19, 3300, 19, 20, 1),
    ("Barringer Crater", 50, 3300, 16, 35, 5),
    ("Tunguska", 60, 3000, 17, 30, 3),
    ("Regional Impact", 150, 3000, 20, 45, 10),
    ("1 km Catastrophic", 1000, 3000, 22, 45, 30),
]

rows = []
for name, D, rho, v, ang, strength in scenarios:
    E_mt = kinetic_energy_mt(D, rho, v)
    h_km = breakup_altitude_km(strength, rho, v)
    crater_km = crater_diameter_km(E_mt if h_km < 1 else 0.1 * E_mt)
    Mw = seismic_magnitude(E_mt * 0.3)
    rows.append({
        "Scenario": name,
        "Diameter (m)": D,
        "Velocity (km/s)": v,
        "Angle (°)": ang,
        "Density (kg/m³)": rho,
        "Strength (MPa)": strength,
        "KE (Mt TNT)": round(E_mt, 2),
        "Breakup Alt (km)": round(h_km, 1),
        "Crater Diam (km)": round(crater_km, 2),
        "Seismic Mw": round(Mw, 2),
        "Event Type": classify_event(h_km)
    })

df = pd.DataFrame(rows)
print(df.to_string(index=False))
