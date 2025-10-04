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
- **NASA SBDB (Small-Body Database Lookup)** — Exercised in `tests/run_data_source_check.py` to cross-check NeoWs values. The lookup endpoint only accepts documented keys (`sstr`, `des`, `spk`) plus hyphenated options like `phys-par` and `full-prec`; the test script now filters requests to that allowlist and flattens the name/value arrays in the response (`orbit.elements`, `phys_par`) so downstream access uses simple dictionaries.
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
| **EUC** | $1.87$ | 20-250 | **Eucrite**: A stony type from the crust of 4 Vesta, similar to volcanic rock (basalt); very strong. | [1], [2] |
| **DIO** | $2.06$ | 20-250 | **Diogenite**: A stony type from deep within the crust of 4 Vesta; also very strong. | [1], [2] |
| **HOW** | $1.89$ | 20-150 | **Howardite**: A stony type from the surface of asteroid 4 Vesta; a mix of Eucrite and Diogenite. | [1], [2] |
| **JAB / IIAB** | $4.46–4.72$ | 100-500 | **Iron Meteorites**: Fragments from the metallic cores of large, molten asteroids; extremely strong. | [1], [2] |

---

### Formulas and Reasoning Used in This Table

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

### Reasoning for Strength Ranges

The reason observational data often gives a minimum value (e.g., "> 20 MPa") is due to a selection bias: these meteoroids are so strong that they typically survive their passage through the atmosphere without breaking up. Scientists can only measure the peak aerodynamic pressure they endured before landing, so their true strength is known to be *at least* that high, but the upper limit remains unmeasured from atmospheric entry alone. The ranges below are established by comparing their composition to analogous materials.

#### **Eucrites & Diogenites (20–250 MPa)**
-   **Observational Lower Limit**: Both Eucrites and Diogenites are strong, crystalline rocks from the asteroid Vesta that have been observed to withstand at least **20 MPa** of pressure during atmospheric entry.
-   **Material Analogue**: Eucrites are basaltic rocks, compositionally similar to volcanic basalt on Earth. The compressive strength of terrestrial basalt typically falls between **100 and 250 MPa**. This provides a reasonable upper bound for their ultimate failure strength.

#### **Howardites (20–150 MPa)**
-   **Composition**: Howardites are a mix of Eucrite and Diogenite fragments that have been fused together into a single rock by impacts on Vesta's surface (a "breccia").
-   **Points of Weakness**: While the individual fragments are very strong, the cement binding them together is a structural weakness. Because of this, the overall strength of a Howardite is likely lower than that of a pure, solid Eucrite or Diogenite.

#### **Iron Meteorites (> 100 MPa, typically 300-500 MPa)**
-   **Observational Lower Limit**: The "> 100 MPa" value is a very conservative lower limit based on the immense pressures they have been observed to survive.
-   **Material Analogue**: Iron meteorites are made of extremely tough iron-nickel alloys (kamacite and taenite). The ultimate tensile and compressive strength of these materials is comparable to high-grade structural steel, often in the range of **300 to 500 MPa**. They are the strongest natural objects known to enter Earth's atmosphere.

---

### Sources

1.  **Mathias, D. L., Wheeler, L. F., & Dotson, J. L. (2017).** *A probabilistic asteroid impact risk model: Assessment of sub-300 m impacts*. Icarus, 289, 106-119. (Density data derived from Table 1 and macroporosity model).
2.  **Popova, O. P., et al. (2011).** *Very low strengths of interplanetary meteoroids and small asteroids*. Meteoritics & Planetary Science, 46(10), 1525-1550. (Strength data sourced from observational analysis).

## Scientific Considerations (Hackathon Guidance)

- **Orbital mechanics** — Model trajectories with Keplerian orbital elements (semi-major axis, eccentricity, inclination, true anomaly) and standard two-body propagation to place Impactor-2025 relative to Earth.
- **Impact energy** — Estimate kinetic energy from mass (derive via size and nominal density ~3000 kg/m³) and approach velocity, then convert to TNT equivalent to communicate scale.
- **Crater scaling** — Apply established scaling relationships (e.g., Collins et al. 2005) to translate impact energy into transient and final crater dimensions.
- **Environmental effects** — Couple impact outcomes to USGS datasets for tsunami inundation (coastal DEM) and seismic shaking proxies, enabling downstream risk visualisation.

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

| White paper parameter | Needed for | Available via NeoWs/SBDB? | Key field(s) | Comments |
| --- | --- | --- | --- | --- |
| Absolute magnitude (*H*) | Diameter inference, size distribution (§2.1) | ✔︎ (NeoWs & SBDB) | `absolute_magnitude_h`, `phys_par.H` | Provided in both feeds |
| Diameter (min/max) | Mass/energy, Monte Carlo sampling (§§2.1–2.3) | ✔︎ (NeoWs & SBDB) | `estimated_diameter.*`, `phys_par.diameter` | NeoWs gives range; SBDB single value/extent |
| Albedo | Converting H to diameter (§2.1) | ✔︎ (SBDB) | `phys_par.albedo` | Covers the NeoWs gap |
| Bulk density | Mass, breakup models (§2.1) | ✔︎ (SBDB when reported) | `phys_par.density` | Fall back to literature/synthetic |
| Material/strength parameters | Fragmentation, breakup (§2.2) | ✘ | — | Use SBDB spectral class as proxy; lab data still needed |
| Entry velocity at Earth | Energy deposition (§2.1) | ✔︎ (NeoWs) | `close_approach_data.relative_velocity` | Units: km/s, km/h, mph |
| Entry angle distribution | Atmospheric entry (§2.1) | ✘ | — | Derived statistically (PAIR Eq. 7) |
| Impact probability / orbit covariance | Scenario weighting (§2.1) | ✘ | — | Requires Sentry/SPICE data |
| Orbital elements (a, e, i, Ω, ω, M) | Trajectory propagation (§2.1) | ✔︎ (NeoWs lookup & SBDB) | `orbital_data.*`, `orbit.*` | Either source works |
| Rotation state, shape | Detailed breakup modeling (§2.2) | ✔︎ (rotation via SBDB) | `phys_par.rot_per` | Shape still sparse |
| Atmospheric profile | Energy deposition (§2.2) | ✘ | — | Use standard atmosphere models |
| Target density, gravity | Crater scaling (§2.3) | ✘ | — | Geological datasets |
| Population distribution | Consequence modeling (§2.5) | ✘ | — | USGS/WorldPop |
| Seismic/tsunami vulnerability curves | Risk metrics (§§2.3–2.4) | ✘ | — | HAZUS/NOAA products |

Only the bolded asteroid-orbit fundamentals (H, diameter, velocity, orbital elements) are covered directly by NeoWs. Environmental, material, and vulnerability parameters must be sourced from SBDB, USGS datasets, or the domain literature acknowledged in the white paper.

## Supplementary Source: NASA Small-Body Database (SBDB)

The SBDB Web Service (`https://ssd-api.jpl.nasa.gov/sbdb.api`) enriches the asteroid record with physical parameters absent in NeoWs. We captured the field mapping in `tests/test_sbdb_data_source.py`, which parses a reference response (Eros 433) and illustrates how to query the live endpoint.

| White paper parameter | Available via SBDB? | SBDB field(s) | Notes |
| --- | --- | --- | --- |
| Absolute magnitude (*H*) | ✔︎ | `phys_par.H` | Redundant with NeoWs |
| Diameter | ✔︎ | `phys_par.diameter` (km), `phys_par.extent` | May differ from NeoWs min/max range |
| Geometric albedo | ✔︎ | `phys_par.albedo` | Covers NeoWs gap |
| Bulk density | ✔︎ (when known) | `phys_par.density` | Not populated for every object |
| Mass / GM | ✔︎ (when known) | `phys_par.mass`, `phys_par.GM` | Useful for momentum considerations |
| Rotation period | ✔︎ | `phys_par.rot_per` | Supports breakup/fragmentation modeling |
| Spectral / taxonomic class | ✔︎ | `phys_par.spec`, `phys_par.spec_B` | Proxy for composition/strength |
| Color indices | ✔︎ | `phys_par.BV`, `phys_par.UB`, etc. | Secondary composition indicators |
| Orbital elements | ✔︎ | `orbit.a`, `orbit.e`, `orbit.i`, `orbit.om`, `orbit.w`, `orbit.ma` | Richer than the basic NeoWs feed |
| Impact probability / covariance | ✘ | — | SBDB does not expose Sentry risk tables |
| Entry angle distribution | ✘ | — | Still needs statistical modeling |
| Environmental layers (DEM, population, tsunami) | ✘ | — | Outside SBDB scope |

Action items:

- Extend a data-ingestion layer to request `phys-par=1&full-prec=1` from SBDB for target objects (e.g., `sstr=99942 Apophis`) and merge the returned albedo, density, rotation, and spectral class into the simulation inputs.
- Fallback to synthetic defaults when SBDB omits a field (e.g., density absent for many small bodies).
- Use the shared designation (`designation` / `object.des`) as the primary key when fusing NeoWs encounter data with SBDB physical parameters.

### Name and identifier crosswalk

The SBDB service uses the MPC/IAU numeric designation (`object.des`), while NeoWs exposes both the designation and a padded `neo_reference_id`. Our exploratory test (`tests/test_sbdb_data_source.py::test_cross_reference_designations`) confirms:

- `SBDB.object.des` == `NeoWs.designation`
- `NeoWs.neo_reference_id` ends with the same numeric designation (e.g., `2000433` → `433`), so string suffix matching or explicit `designation` should be used when joining the two sources.

Use the shared designation as the primary key when fusing NeoWs encounter data with SBDB physical parameters.

For a quick terminal report of which PAIR parameters map to NeoWs versus SBDB, run `python3 tests/print_data_source_matrix.py`; this prints the same availability matrix captured above and can optionally fetch live SBDB data when `LIVE=1` is set (override the target with `SBDB_ID=Apophis` or a preferred designation).

To exercise the full integration pipeline, `python3 tests/run_data_source_check.py` attempts to pull the white-paper parameters from both APIs (it auto-loads `.env` for keys) and prints the retrieved values side by side, aborting if either service is unreachable.
