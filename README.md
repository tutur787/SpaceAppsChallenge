# SpaceAppsChallenge

Meteor Madness is a Streamlit MVP that demonstrates how live NASA NeoWs data can seed an educational asteroid-impact simulator.

## Requirements

- Python 3.9+
- `streamlit`, `numpy`, `pandas`, `requests`, `pydeck`
- Optional: `NASA_API_KEY` environment variable (falls back to `DEMO_KEY` for low-rate access)

To install dependencies once:

```bash
python3 -m pip install --upgrade streamlit numpy pandas requests pydeck
```

## Running the app

```bash
export NASA_API_KEY=<your_api_key>  # optional
python3 -m streamlit run main.py
```

On startup `main.py` fetches the current-day NeoWs feed. The first hazardous object (or the first object if no hazardous entries exist) becomes the **baseline** that hydrates the diameter and velocity sliders. If the API call fails, the app raises a warning and reverts to the documented whitepaper defaults.

## Telemetry & provenance

- The sidebar now lists dataset provenance (NeoWs feed status), slider validation results, and which sliders received live values.
- For headless QA, set `METEOR_MADNESS_HEADLESS_TELEMETRY=1` before launching the app. The same telemetry summary is printed to stdout so you can capture it in logs/CI.

## Tests

Run the regression suite, which includes slider spec guards and telemetry helpers:

```bash
python3 -m pytest
```

The NeoWs integration check remains available via `python3 -m unittest tests.test_nasa_neo_live` when a live API key is configured.
