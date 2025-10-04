# Data Sources

This file links the data requirements outlined in Mathias et al. (2017) — the white paper that inspired the PAIR framework — to the datasets the MVP currently uses or stubs. Columns labelled *NeoWs?* indicate whether the NASA Near-Earth Object Web Service (NeoWs) already exposes the parameter.

## Asteroid and Entry Parameters

| Model component (PAIR §) | Required fields | Preferred source | NeoWs? | Status in MVP |
| --- | --- | --- | --- | --- |
| Size distribution (§2.1) | Diameter (m), absolute magnitude *H*, albedo | NASA NeoWs `estimated_diameter`, `absolute_magnitude_h` | ✔︎ (`estimated_diameter`, `absolute_magnitude_h`) | Diameter min/max pulled; *H* not yet stored |
| Density / composition (§2.1) | Bulk density, material class | Taxonomy catalogues (e.g., SMASS), literature | ✘ | Synthetic presets in `main.py` |
| Strength scaling (§2.2) | Bulk compressive strength, fragmentation exponent α | Laboratory / literature (white paper Table 1) | ✘ | Slider with synthetic defaults |
| Impact velocity (§2.1) | Relative velocity at encounter (km/s) | NeoWs `close_approach_data.relative_velocity` | ✔︎ | Used in Explore tab |
| Entry angle distribution (§2.1) | Entry angle θ | Derived statistically (PAIR Eq. 7) | ✘ | Monte Carlo sampler approximates distribution |
| Orbit / ephemeris | Semi-major axis, eccentricity, inclination | NeoWs `orbital_data` | ✔︎ | Not yet fetched; needed for orbit visualizer |

## Atmospheric Interaction & Energy Deposition

| Model component | Required fields | Preferred source | NeoWs? | Status in MVP |
| --- | --- | --- | --- | --- |
| Atmospheric density profile (§2.2) | ρ_air vs altitude | US Standard Atmosphere / NRLMSISE-00 | ✘ | Simplified exponential constants baked in |
| Breakup criterion (§2.2) | Aerodynamic pressure threshold | White paper Table 1 parameters | ✘ | Uses slider strength and exponential atmosphere |
| Energy deposition history (§2.2) | dE/dh, airburst altitude | Derived from entry solver | ✘ | Simplified ground-coupling fraction |

## Surface Effects & Consequences

| Model component | Required fields | Preferred source | NeoWs? | Status in MVP |
| --- | --- | --- | --- | --- |
| Blast overpressure (§2.3) | Threshold curves (1, 4, 12 psi) | White paper Eq. 5, Glasstone & Dolan | ✘ | Implemented analytically |
| Thermal radiation (§2.3) | Luminous efficiency η, heat flux thresholds | White paper Eq. 6 | ✘ | Not yet modelled |
| Seismic coupling (§2.3) | Coupling fraction κ, magnitude relation | White paper Eq. 8 | ✘ | Uses constant κ = 3e-4 |
| Crater scaling (§2.3) | π-group constants, transition rules | Collins et al. 2005 / white paper refs | ✘ | Simplified transient-to-final crater logic |
| Tsunami / bathymetry (§2.4) | Bathymetry, coastal DEM | USGS National Map DEM, NOAA SRTM30+ | ✘ | Not yet integrated |
| Population exposure (§2.5) | Population raster, vulnerability curves | USGS/HAZUS, WorldPop, LandScan | ✘ | Synthetic densities & casualty fractions |

## External Datasets Currently Wired

- **NASA NeoWs (Near-Earth Object Web Service)** — Queried through `fetch_today_neos()`. Currently harvests asteroid name, diameter range, hazardous flag, relative velocity, and miss distance. Additional fields such as `absolute_magnitude_h` and `orbital_data` are available via the same endpoint but are not parsed yet; add them when we build orbit visualisations or size-from-magnitude conversions.
- **Basemaps** — Mapbox Dark (with `MAPBOX_API_KEY`) or OpenStreetMap tiles for context maps.

## Synthetic Placeholders

- **Concentric population densities** — Severe: 5,000 people/km²; Moderate: 2,000; Light: 800. See `SYNTHETIC_POP_DENSITY` in `main.py`.
- **Casualty fractions** — Severe 35%, Moderate 10%, Light 2% (`SYNTHETIC_CASUALTY_RATE`).
- **Atmosphere** — Exponential scale height (7.16 km) and sea-level density constants approximate the profile referenced in Mathias et al.
- **Monte Carlo priors** — Log-normal spreads around user-selected size, density, strength, and the canonical entry-angle distribution (PAIR Eq. 7) encoded in `run_pair_simulation()`.

Replace the synthetic elements with the cited datasets (USGS DEMs, WorldPop/HAZUS rasters, SMASS taxonomy tables, etc.) before using the tool for decision support.
