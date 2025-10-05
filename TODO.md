# TODO

- [x] Align damage/population model with PAIR methodology (Mathias et al. 2017) - COMPLETED
    1. ✅ Renamed expander from "Synthetic exposure breakdown" to "PAIR damage assessment (Mathias et al. 2017)" (`main.py:549`).
    2. ✅ Updated expander to highlight 4-psi blast overpressure as primary damage threshold (white paper Section 2.3, line 236).
    3. ✅ Replaced hardcoded casualty estimates with PAIR "affected population" metric (everyone within 4-psi circle).
    4. ✅ Main metrics now show: "PAIR affected population (≤4-psi)" instead of "Estimated casualties" (`main.py:539`).
    5. ✅ Added documentation that `SYNTHETIC_CASUALTY_RATE` is NOT from PAIR; kept for legacy purposes only (`main.py:58-66`).
    6. ✅ Fixed calculation: 4-psi is a full circle (area = π × r²), not a ring (annulus) (`main.py:531, 559`).
    7. ✅ Added UI warning: "Affected ≠ casualties; depends on building codes, warning time, shelter availability."
    8. ✅ Kept three blast rings (12/4/1 psi) for visualization but emphasized 4-psi as PAIR standard.
    9. 🔲 (Future) Implement Collins et al. (2005) thermal radiation damage radius (white paper Eq. 6, Section 2.4).
    10. 🔲 (Future) Replace synthetic population densities with SEDAC gridded census data (2.5-arc-minute resolution, white paper Section 2.5).

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
