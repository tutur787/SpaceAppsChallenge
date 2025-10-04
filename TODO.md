# TODO

- [ ] Replace synthetic asteroid material presets with SBDB-derived defaults.
    1. Update the NeoWs ingestion (`main.py::fetch_today_neos` or a new `src/api` helper) to join each object with `https://ssd-api.jpl.nasa.gov/sbdb.api?sstr=<designation>&phys-par=1&full-prec=1` so the live payload contains `phys_par.spec`, `phys_par.albedo`, `phys_par.density`, and `phys_par.rot_per`.
    2. Drop those fields into a cache (e.g., `st.session_state['sbdb_cache']`) keyed by designation so Streamlit widgets can pre-populate without a second network call.
    3. Build a lookup table for taxonomy -> {density, strength} using Mathias et al. (2017) Table 1 (e.g., C-complex ~ 1400 kg/m^3, S-complex ~ 3300 kg/m^3, X-complex ~ 4200 kg/m^3). Store it next to `MATERIAL_PRESETS`.
    4. When `phys_par.density` is missing, compute density/strength from the taxonomy lookup and overwrite the default slider values before rendering the UI (`st.session_state['density_<material>'] = ...`).
    5. Log the provenance ("SBDB phys_par" vs "taxonomy fallback") in the sidebar so it is obvious whether the value is measured or inferred.

- [ ] Implement crater scaling with Collins et al. (2005) pi-group equations driven by real target properties.
    1. Create `compute_transient_crater(diameter_m, velocity_km_s, density_kg_m3, impact_angle_deg, target_density_kg_m3)` in `main.py` (or `src/physics.py`) encoding the gravity- and strength-dominated regimes from Collins Table 1.
    2. Replace `final_crater_diameter_m` with the new helper and preserve the transient->final correction (simple vs complex) cited in the white paper.
    3. Source `target_density_kg_m3` from terrain: call the USGS National Map Elevation Point Query Service for the selected lat/lon and branch to seawater (~1025 kg/m^3) when elevation < 0, otherwise crustal rock (~2500 kg/m^3). Cache the DEM response alongside the SBDB cache.
    4. Add unit tests in `tests/test_crater_scaling.py` comparing Chelyabinsk-size and Chicxulub-size benchmarks against published results (+/-20% tolerance).

- [ ] Model environmental effects with USGS data layers so tsunami and seismic outputs stop using constants.
    1. Wire `fetch_neic_events(lat, lon, radius_km, start_date, end_date)` (USGS NEIC API) and store magnitude-distance pairs. Use them to calibrate `DEFAULT_SEISMIC_COUPLING` by fitting Eq. 8 from Mathias et al. (2017).
    2. For tsunami depth, request a 1 arc-second DEM tile around the impact point (USGS National Map WMS). Sample bathymetry/elevation along radial transects (e.g., every 10 deg out to 200 km) and attenuate wave height using NOAA run-up curves.
    3. Surface those layers in the Explore tab: add a `pydeck.Layer` heatmap for peak ground acceleration and a filled contour for expected inundation depth.
    4. Document the API endpoints, query parameters, and update cadence in `DATA_SOURCES.md` so downstream contributors know how to refresh credentials.

- [ ] Propagate orbital elements for Impactor-2025 so trajectory visuals are physics-backed.
    1. After merging NeoWs/SBDB records, retain `orbit.{epoch_tdb,a,e,i,om,w,ma}` in the cached payload.
    2. Implement `propagate_orbit(orbit_elements, target_epoch)` that converts mean anomaly -> eccentric anomaly -> true anomaly, returning Earth-centered Cartesian coordinates (two-body solution is sufficient for the MVP).
    3. Render the propagated arc with `pydeck.LineLayer`/`ScatterplotLayer`, anchoring Earth at the origin and reusing the propagated state to drive inbound bearing in the Defend tab.
    4. Add `tests/test_orbit_propagation.py` that propagates Apophis to 2029-04-13 and asserts the miss distance matches JPL Horizons to within +/-1%.

- [ ] Update docs so contributors know how to run the new data fetches.
    1. Extend `README.md` with a "Data prerequisites" section covering NASA API key, SBDB usage, USGS tokens, and caching strategy.
    2. Note the new optional `.env` keys (`USGS_API_KEY`, `SBDB_CACHE_PATH`) and reference them in `tests/run_data_source_check.py`.
    3. Embed quickstart snippets for running the crater/tsunami unit tests and the updated data-source check.

- [ ] Force MVP sliders to mirror whitepaper parameters exactly.
    1. Build a slider-to-parameter mapping table from `DATA_SOURCES.md` that enforces identical labels, units, bounds, and default step sizes.
        - 2025-10-04: `src/config/slider_specs.py` created and the Streamlit UI now consumes the canonical specs.
    2. Add validation in the UI layer that raises/logs when a slider diverges from the whitepaper spec so regressions are obvious during development.
        - 2025-10-04: `_validate_slider_args()` in `main.py` now fails fast (and surfaces `st.error`) whenever slider kwargs drift from the spec table.
    3. Document any intentional deviations (if required) inline in `DATA_SOURCES.md` so reviewers know the rationale.

- [ ] Seed slider defaults with live API values instead of hard-coded numbers.
    1. Fetch baseline values from NeoWs/SBDB before rendering Streamlit widgets and hydrate `st.session_state` accordingly.
        - 2025-10-04: Explore tab seeds defaults from today's NeoWs feed with fallback warnings when live data is unavailable.
    2. Provide a fallback path (with an explicit warning) if the API call fails so we never silently fall back to synthetic presets.

- [ ] Detect and flag any synthetic data in the MVP pipeline.
    1. Tag every dataset record with a provenance flag (`live`, `cached`, `synthetic`) during ingestion.
    2. Surface the flag in both the logging/CLI printout and the frontend sidebar so users know when data is fabricated.
        - 2025-10-05: Sidebar + headless telemetry scaffolding landed (`src.telemetry`) with provenance labels (`live`, `material_preset`, `whitepaper`, `spec_floor`). Extend to enforce synthetic detection in Phase 2.
    3. Extend telemetry to differentiate cached vs synthetic feeds and fail the headless report when synthetic data appears unexpectedly (Phase 2 follow-up).

- [ ] Audit dropdown behaviour for whitepaper and synthetic views.
    1. Write a smoke test that asserts the "White paper inputs" dropdown emits the expected options instead of rendering blank.
    2. Decide whether the "Synthetic exposure breakdown" dropdown should exist; if removal is chosen, excise the component and associated state.

- [ ] Replace the "Today's NEO data" table with a selector-driven workflow.
    1. Render a dropdown (or searchable select) populated with the day's NeoWs objects; selecting one should hydrate the detail panels.
    2. Preserve access to the full table via a secondary "View details" action so power users can still inspect the dataset.

- [ ] Print a full telemetry dump when running `main.py` non-interactively.
    1. Emit a terminal report summarising slider bindings, data provenance flags, and dropdown option counts so QA can run in headless mode.
    2. Fail fast (non-zero exit) when the report detects missing bindings, blank dropdowns, or synthetic data where it should not appear.
        - 2025-10-05: Added environment-controlled headless summary via `METEOR_MADNESS_HEADLESS_TELEMETRY`; expand with exit codes and dropdown metrics in Phase 4.
