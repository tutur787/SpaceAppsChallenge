# Continuation Prompt for Next LLM

## Repository Snapshot (2025-10-05)
- Branch/worktree status: slider fidelity (Phase 1) hardening in progress; telemetry scaffolding freshly added.
- Key modules:
  - `main.py` — Streamlit UI; now emits provenance/validation telemetry via `render_sidebar_telemetry()` and optional headless logging gated by `METEOR_MADNESS_HEADLESS_TELEMETRY`.
  - `src/telemetry.py` — new helpers for provenance tags, slider/dataset status tracking, and formatted headless reports.
  - `tests/test_telemetry.py` — unit coverage for telemetry helpers.
- Docs updated: `README.md`, `DATA_SOURCES.md`, `TODO.md`, `ROADMAP.md`.

## What’s Done
1. Slider defaults now record provenance (`live`, `material_preset`, `whitepaper`, `spec_floor`, etc.) and surface status in the sidebar/headless output.
2. Dataset ingestion (`load_today_neo_feed`, `seed_live_slider_defaults`) logs NeoWs feed health and hydrates telemetry state.
3. Headless telemetry toggle via `METEOR_MADNESS_HEADLESS_TELEMETRY` prints the consolidated report for CI/QA.
4. Documentation refreshed to explain NeoWs hydration workflow and telemetry plumbing.
5. Tests: `tests/test_telemetry.py` green; entire suite runs but fails on `tests/test_nasa_neo_live` due to sandbox DNS blocks when resolving `api.nasa.gov` (see latest pytest output).

## Outstanding / Next Steps
1. **Phase 2 provenance expansion**
   - Differentiate cached vs synthetic sources (extend `ProvenanceTag` usage).
   - Enforce warnings/errors when synthetic data appears in live mode.
   - Add targeted tests (mocked ingestion) to validate enforcement logic.
2. **Dropdown hygiene (Phase 3 prep)**
   - Review "Synthetic exposure breakdown" expander and the whitepaper dropdown workflow.
   - Plan telemetry hooks to count dropdown options for future headless checks.
3. **Headless exit codes (Phase 4)**
   - Decide on CLI entry point (`python main.py --headless-report` or similar) and integrate telemetry summary with non-zero exits when mismatches occur.
4. **Fix live NeoWs tests in sandbox**
   - Requires network allowance or mocking; currently blocked by DNS errors.

## How to Resume
1. Install deps if needed: `python3 -m pip install --upgrade streamlit numpy pandas requests pydeck pytest`.
2. Optional: export `NASA_API_KEY` to hydrate live defaults (Streamlit runtime); tests fall back to `DEMO_KEY`.
3. Run app: `python3 -m streamlit run main.py` (sidebar shows telemetry summary).
4. Run tests: `python3 -m pytest` (expect NeoWs live tests to fail without network — safe to skip with `-k "not nasa_neo_live"`).
5. When extending telemetry, update `DATA_SOURCES.md` and `TODO.md` with any intentional spec deviations or provenance changes.

## Handoff Notes
- Maintain consistency between `SLIDER_SPECS` and UI helper usage; `_validate_slider_args()` now logs outcomes and should continue to be the single source of spec validation.
- Telemetry report stored in `st.session_state["telemetry_report"]`; reuse `ensure_report()` when adding new instrumentation.
- When adding provenance for new datasets (e.g., USGS), call `record_dataset_provenance()` with clear details so the sidebar/headless views stay concise.

