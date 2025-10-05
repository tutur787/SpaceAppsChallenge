"""Canonical slider specifications derived from the PAIR white paper.

This module captures the authoritative parameter bounds and defaults cited in
`DATA_SOURCES.md` and `White_paper.md` so Streamlit widgets can stay aligned
with the literature. Phase 1 will replace the inline widget arguments in
`main.py` with these specs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional


@dataclass(frozen=True)
class SliderSpec:
    """Metadata describing how a UI slider should behave."""

    key: str
    label: str
    unit: str
    min_value: float
    max_value: float
    step: float
    whitepaper_default: Optional[float]
    whitepaper_source: str
    description: str = ""
    requires_live_default: bool = True
    ui_context: str = "explore"


SLIDER_SPECS: Dict[str, SliderSpec] = {
    # Mathias et al. (2017) restricts the ensemble study to impactors <= 300 m
    # and highlights 65 m as the risk threshold (§3.1, Fig. 12). Bounds taken
    # from DATA_SOURCES.md “White Paper Parameters vs. NeoWs Coverage”.
    "diameter_m": SliderSpec(
        key="diameter_m",
        label="Asteroid diameter (m)",
        unit="m",
        min_value=20,
        max_value=300,
        step=5,
        whitepaper_default=65,
        whitepaper_source="White_paper.md:329-339",
        description=(
            "PAIR samples diameters derived from H=20–30 and albedo distributions;"
            " 300 m is the upper bound of the baseline study."
        ),
    ),
    # Impact velocities from the Greenstreet et al. (2012) distribution (White_paper.md:415-420)
    # cluster around 11–35 km/s for Earth encounters. We clamp to that physically
    # representative range.
    "velocity_km_s": SliderSpec(
        key="velocity_km_s",
        label="Velocity at impact (km/s)",
        unit="km/s",
        min_value=11.0,
        max_value=35.0,
        step=0.5,
        whitepaper_default=19.0,
        whitepaper_source="White_paper.md:415-420",
        description=(
            "Greenstreet et al. (2012) NEO ensemble feeding PAIR spans roughly 11–35 km/s;"
            " 19 km/s is near the median encounter speed."
        ),
    ),
    # Entry angles are sampled from a cosine distribution on 0°–90° with a mode at 45° (Shoemaker 1961).
    "angle_deg": SliderSpec(
        key="angle_deg",
        label="Impact angle (°)",
        unit="deg",
        min_value=0,
        max_value=90,
        step=1,
        whitepaper_default=45,
        whitepaper_source="White_paper.md:417-423",
        description=(
            "White paper samples entry angles from a cosine distribution favouring 45°."
        ),
        requires_live_default=False,
        ui_context="monte_carlo",
    ),
    # Density limits from meteorite distribution with macroporosity clipping to 1.1–7.5 g/cm³ (White_paper.md:380-402).
    "bulk_density": SliderSpec(
        key="bulk_density",
        label="Bulk density (kg/m³)",
        unit="kg/m³",
        min_value=1100,
        max_value=7500,
        step=50,
        whitepaper_default=2500,
        whitepaper_source="White_paper.md:380-402",
        description=(
            "Meteorite-derived density distribution trimmed to macroporosity bounds;"
            " 2.5 g/cm³ is the PAIR baseline for stony bodies."
        ),
        requires_live_default=False,
        ui_context="defend",
    ),
    # Aerodynamic breakup strength sampled logarithmically from 0.1–10 MPa (White_paper.md:404-411).
    "strength_mpa": SliderSpec(
        key="strength_mpa",
        label="Bulk compressive strength (MPa)",
        unit="MPa",
        min_value=0.1,
        max_value=10.0,
        step=0.1,
        whitepaper_default=1.0,
        whitepaper_source="White_paper.md:404-411",
        description=(
            "Popova et al. (2011) bounds used in PAIR for initial breakup strength;"
            " logarithmic sampling spans 0.1–10 MPa."
        ),
        requires_live_default=False,
        ui_context="defend",
    ),
    # Monte Carlo sampling size: MVP uses 100–1000 for interactivity; document that this is a downscaled
    # subset of the 30M draws in the white paper.
    "pair_samples": SliderSpec(
        key="pair_samples",
        label="Number of Monte Carlo samples",
        unit="count",
        min_value=100,
        max_value=1000,
        step=50,
        whitepaper_default=300,
        whitepaper_source="White_paper.md:432-448",
        description=(
            "UI exposes a reduced sample count for responsiveness while referencing the"
            " 30M-case study in Mathias et al."
        ),
        requires_live_default=False,
        ui_context="defend",
    ),
    # Deflection controls: not from the PAIR paper, but we retain bounds consistent with educational delta-v studies.
    "delta_v_mm_s": SliderSpec(
        key="delta_v_mm_s",
        label="Δv (mm/s)",
        unit="mm/s",
        min_value=0.0,
        max_value=5.0,
        step=0.1,
        whitepaper_default=0.5,
        whitepaper_source="DefendTab",
        description="Educational deflection slider; bounds informed by kinetic-impactor examples.",
        requires_live_default=False,
    ),
    "lead_days": SliderSpec(
        key="lead_days",
        label="Lead time (days)",
        unit="days",
        min_value=0,
        max_value=3650,
        step=30,
        whitepaper_default=365,
        whitepaper_source="DefendTab",
        description="Lead time range covers up to 10 years for what-if planning.",
        requires_live_default=False,
    ),
    "inbound_bearing": SliderSpec(
        key="inbound_bearing",
        label="Inbound bearing (°)",
        unit="deg",
        min_value=0.0,
        max_value=359.0,
        step=1.0,
        whitepaper_default=90.0,
        whitepaper_source="DefendTab",
        description="Bearing is a full azimuth sweep for the simplified deflection map.",
        requires_live_default=False,
    ),
}


def get_slider_spec(key: str) -> SliderSpec:
    """Return the slider spec for `key`, raising a helpful error if missing."""

    try:
        return SLIDER_SPECS[key]
    except KeyError as exc:  # pragma: no cover - defensive guard
        known = ", ".join(sorted(SLIDER_SPECS))
        raise KeyError(f"Unknown slider spec '{key}'. Known keys: {known}") from exc
