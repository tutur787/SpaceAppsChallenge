# SBDB Lookup Failure – Handoff Notes

## Overview
The `tests/run_data_source_check.py` script cross-references NASA NeoWs and JPL SBDB data for an asteroid (default: Apophis). NeoWs fetches succeed, but the SBDB request now returns HTTP 400.

## Current Behavior
Running `python3 tests/run_data_source_check.py` produces:

```
ERROR: SBDB request failed after attempts: 400 Client Error: Bad Request ...
{"code":"400","message":"one or more query parameter was not recognized"}
```

The script tries several SBDB query variants (designation number, SPK id, full name, tail token) but every response comes back with the same error.

## Environment Notes
- Python 3.9 on macOS, `requests`/`urllib3` available.
- Outbound HTTP to `api.nasa.gov` and `ssd-api.jpl.nasa.gov` is permitted when the user runs the script (the local sandbox environment blocks it, so remote runs are necessary for verification).
- `NASA_API_KEY` is provided via `API.env`.

## Work So Far
- Updated NeoWs default ID to the current reference id for Apophis (`3542519`).
- Added richer SBDB logging; we now capture the full JSON error payload for 4xx responses.
- Adjusted SBDB queries to use `phys=1` instead of the old `phys-par=1` flag, per the latest API docs (see `tests/run_data_source_check.py:97-117`).
- Documentation and helper utilities (`tests/test_sbdb_data_source.py`, `tests/print_data_source_matrix.py`, `DATA_SOURCES.md`) were updated to match the new `phys` flag.

Despite the switch to `phys=1`, SBDB still replies with `"one or more query parameter was not recognized"`. That means at least one of the parameters we are sending is incorrect.

## Hypotheses / Next Steps
1. Re-read the current SBDB API documentation and confirm the correct query-string keys (the API might have renamed `orb` and/or `phys`).
2. Test a minimal request with only the required identifier parameter (e.g., `https://ssd-api.jpl.nasa.gov/sbdb.api?s=Apophis`) to confirm baseline connectivity, then add optional flags one at a time.
3. If the API now expects `full-phy` or other flags, adjust `build_sbdb_attempts` accordingly.
4. Ensure the default identifier matches SBDB’s expectations (for Apophis, they might prefer `s=Apophis` or `spk=2000423`).

## Deliverable Expectations
- Script should exit successfully, printing the comparison table using live SBDB data (crater/physical parameters populated).
- Tests that rely on the SBDB helper (`tests/test_sbdb_data_source.py`) should pass with the updated parameter names.

Please continue the investigation to identify the new parameter key(s) the SBDB API accepts and update the client accordingly.
