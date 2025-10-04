# Data Sources

## Live / External

- **NASA NeoWs (Near Earth Object Web Service)** — Queried through `fetch_today_neos()` using the NASA demo key or an environment variable `NASA_API_KEY`. This supplies real asteroid diameter and velocity estimates when network access is permitted.
- **Basemaps** — Mapbox Dark style is used if the `MAPBOX_API_KEY` is present; otherwise an OpenStreetMap tile layer is requested. Both are external providers and require network access to render maps.

## Synthetic Placeholders

- **Concentric population densities** — When no geospatial population raster is available, the app assumes the following uniform densities (people/km²) for the concentric damage zones:
  - Severe (≥12 psi): 5,000
  - Moderate (≈4 psi): 2,000
  - Light (≈1 psi): 800
  These values live in `SYNTHETIC_POP_DENSITY` inside `main.py` and stand in for USGS/WorldPop rasters. They represent dense urban, suburban, and rural gradients respectively.
- **Casualty fractions** — To translate exposed population into casualties, the MVP applies fixed ratios: 35% (severe), 10% (moderate), 2% (light). These are intentionally conservative placeholders derived from historical blast injury studies; replace them with location-specific vulnerability curves when real demographic data are wired in.
- **Monte Carlo priors** — The probabilistic simulator samples log-normal spreads around user-selected diameter, density, and strength because detailed uncertainty distributions are not bundled with the MVP. Documented in code comments referencing Mathias et al. (2017).

All synthetic constants are documented inside the repository and should be replaced with authoritative datasets (e.g., USGS National Map elevation, population grids, or HAZUS vulnerability relationships) before operational use.

## Notes

- Atmospheric breakup heights, blast-radius scaling, and crater-sizing follow published relationships summarized in Mathias et al. (2017) and Collins et al. (2005). While the formulas are physics-based, they are still approximations; validate against high-fidelity tools when integrating mission data.
- Deflection outcomes reuse the same synthetic exposure model because site-specific terrain and population rasters are not yet integrated.
