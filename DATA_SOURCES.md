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

## Asteroid Classification, Density, and Strength (Updated)

This revised table provides the estimated **bulk density** of asteroids by applying a 34% porosity to meteorite samples. The aerodynamic strength ranges are updated based on observational data for each specific classification.

| Classification | Estimated Bulk Density (g/cm³) | Aerodynamic Strength Range (MPa) | Meaning | Source (Density, Strength) |
| :--- | :--- | :--- | :--- | :--- |
| **CM** | $1.50$ | 0.01–0.1 | **Carbonaceous (Mighei-like)**: A carbon-rich type containing water and organic compounds; very fragile. | [1], [2] |
| **LL** | $2.11$ | 0.1–1 | **Low Iron, Low Metal**: A stony type with low total iron and very little in metallic form. | [1], [2] |
| **L** | $2.18$ | 1–4 | **Low Iron**: A common stony type with a lower iron content than H-types. | [1], [2] |
| **H** | $2.23$ | 4–20 | **High Iron**: A common stony type with a high amount of metallic iron; relatively strong. | [1], [2] |
| **EUC** | $1.87$ | > 20 | **Eucrite**: A stony type from the crust of 4 Vesta, similar to volcanic rock (basalt); very strong. | [1], [2] |
| **DIO** | $2.06$ | > 20 | **Diogenite**: A stony type from deep within the crust of 4 Vesta; also very strong. | [1], [2] |
| **HOW** | $1.89$ | > 20 | **Howardite**: A stony type from the surface of asteroid 4 Vesta; a mix of Eucrite and Diogenite. | [1], [2] |
| **JAB / IIAB** | $4.46–4.72$ | > 100 | **Iron Meteorites**: Fragments from the metallic cores of large, molten asteroids; extremely strong. | [1], [2] |

---

### Formulas Used in This Table

#### Calculating Bulk Density

The bulk density of an asteroid is estimated by adjusting the solid (grain) density of its corresponding meteorite type for macroporosity (empty space).

The formula is:
$$\rho_{bulk} = \rho_{grain} \times (1 - \phi)$$

-   **$\rho_{bulk}$**: The bulk density of the entire asteroid.
-   **$\rho_{grain}$**: The grain density of the solid meteorite material.
-   **$\phi$**: The macroporosity, or the fraction of empty space (a mean value of 0.34, or 34%, is used here).

**Example (H-type asteroid):**
-   Grain density = $3.38 \text{ g/cm}^3$
-   Calculation: $3.38 \times (1 - 0.34) = 3.38 \times 0.66 = 2.23 \text{ g/cm}^3$.

#### Defining Aerodynamic Strength

Aerodynamic strength is the maximum **ram pressure** an object can withstand before breaking up during atmospheric entry. The values in the table are the measured strengths at which this breakup occurs.

The formula for ram pressure is:
$$P_{ram} = \rho_{air} v^2$$

-   **$P_{ram}$**: The ram pressure, the force exerted on the asteroid by the atmosphere.
-   **$\rho_{air}$**: The density of the air at the asteroid's altitude.
-   **$v$**: The velocity of the asteroid.

An asteroid breaks apart when **$P_{ram}$ > Aerodynamic Strength**.

---

### Sources

1.  **Mathias, D. L., Wheeler, L. F., & Dotson, J. L. (2017).** *A probabilistic asteroid impact risk model: Assessment of sub-300 m impacts*. Icarus, 289, 106-119. (Density data derived from Table 1 and macroporosity model).
2.  **Popova, O. P., et al. (2011).** *Very low strengths of interplanetary meteoroids and small asteroids*. Meteoritics & Planetary Science, 46(10), 1525-1550. (Strength data sourced from observational analysis).

## NeoWs Field Coverage

The NeoWs feed returns a rich set of attributes for each close-approach record. The table below marks which columns the API exposes and whether the current MVP ingests them.

| NeoWs field | Description | Present in API response | Used in MVP? | Notes / next steps |
| --- | --- | --- | --- | --- |
| `id`, `neo_reference_id` | Internal identifiers | ✔︎ | ✘ | Store when we add caching or detail drill-down |
| `name` | Asteroid designation | ✔︎ | ✔︎ | Displayed in optional NEO table |
| `absolute_magnitude_h` | H-magnitude (brightness proxy for size) | ✔︎ | ✘ | Needed to derive diameter when only H is known |
| `estimated_diameter_min_km`, `estimated_diameter_max_km` (plus miles/ft variants) | Diameter estimates per unit | ✔︎ | ✘ | We currently read only the meter values; convert on demand |
| `estimated_diameter_min_m`, `estimated_diameter_max_m` | Diameter estimates (meters) | ✔︎ | ✔︎ | Used for context list |
| `is_potentially_hazardous_asteroid` | Hazard flag | ✔︎ | ✔︎ | Surface-level risk indicator |
| `close_approach_date`, `close_approach_date_full`, `epoch_date_close_approach` | Timing of approach | ✔︎ | ✘ | Add to timeline visualisation |
| `relative_velocity` (`kilometers_per_second`, `kilometers_per_hour`, `miles_per_hour`) | Encounter speed | ✔︎ | ✔︎ (`kilometers_per_second`) | Other units available if needed |
| `miss_distance` (`kilometers`, `lunar`, `astronomical`, `miles`) | Miss distance in multiple units | ✔︎ | ✔︎ (`kilometers`) | Lunar/AU handy for outreach graphics |
| `orbiting_body` | Body the asteroid is passing | ✔︎ | ✘ | Useful for multi-body context |
| `links` (`self`, etc.) | Hyperlinks to detail endpoints | ✔︎ | ✘ | Use for deep-dive panels |
| `orbital_data` (present in lookup endpoint) | Keplerian elements, ephemeris data | ✔︎ (separate endpoint) | ✘ | Required for orbit propagation view |

To expand the MVP, extend `fetch_today_neos()` to persist the unused columns (e.g., `absolute_magnitude_h`, `close_approach_date`, `orbiting_body`) and surface them in the UI or downstream calculations (size-from-H conversions, encounter timelines, orbital tracks).

## White Paper Parameters vs. NeoWs Coverage

The PAIR framework (Mathias et al. 2017) calls for the following inputs. This table flags whether each parameter is obtainable directly from the NeoWs API; downstream implementation is tracked elsewhere.

| White paper parameter | Needed for | Available via NeoWs? | NeoWs field(s) if available | Comments |
| --- | --- | --- | --- | --- |
| Absolute magnitude (*H*) | Diameter inference, size distribution (§2.1) | ✔︎ | `absolute_magnitude_h` | Provided in each NEO object |
| Diameter (min/max) | Mass/energy, Monte Carlo sampling (§§2.1–2.3) | ✔︎ | `estimated_diameter` (meters/kilometers) | Multiple unit variants |
| Albedo | Converting H to diameter (§2.1) | ✘ | — | Needs SMASS/WISE catalogues |
| Bulk density | Mass, breakup models (§2.1) | ✘ | — | Literature/synthetic presets |
| Material/strength parameters | Fragmentation, breakup (§2.2) | ✘ | — | Requires lab data; not in NeoWs |
| Entry velocity at Earth | Energy deposition (§2.1) | ✔︎ | `close_approach_data.relative_velocity` | Units: km/s, km/h, mph |
| Entry angle distribution | Atmospheric entry (§2.1) | ✘ | — | Derived statistically (PAIR Eq. 7) |
| Impact probability / orbit covariance | Scenario weighting (§2.1) | ✘ | — | Requires Sentry/SPICE data |
| Orbital elements (a, e, i, Ω, ω, M) | Trajectory propagation (§2.1) | ✔︎ (via lookup endpoint) | `orbital_data` | Need separate NeoWs lookup call |
| Rotation state, shape | Detailed breakup modeling (§2.2) | ✘ | — | Often unknown; not in NeoWs |
| Atmospheric profile | Energy deposition (§2.2) | ✘ | — | Use standard atmosphere models |
| Target density, gravity | Crater scaling (§2.3) | ✘ | — | Geological datasets |
| Population distribution | Consequence modeling (§2.5) | ✘ | — | USGS/WorldPop |
| Seismic/tsunami vulnerability curves | Risk metrics (§§2.3–2.4) | ✘ | — | HAZUS/NOAA products |

Only the bolded asteroid-orbit fundamentals (H, diameter, velocity, orbital elements) are covered by NeoWs. Environmental, material, and vulnerability parameters must be sourced from USGS datasets or domain literature acknowledged in the white paper.
