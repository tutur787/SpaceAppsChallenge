# SpaceAppsChallenge
# How to run:
  - Install (once): python3 -m pip install --upgrade streamlit numpy pandas requests
  pydeck
  - Optional for the NASA feed: export NASA_API_KEY=... (or edit API.env / .env)
  - Launch the app: python3 -m streamlit run main.py
  - To sanity-check without the UI, run a syntax pass: python3 -c "import ast;
  ast.parse(open('main.py').read())"

  If you prefer isolation, create and activate a virtualenv first (python3 -m venv .venv
  && source .venv/bin/activate) and run the same install/start commands inside it.


## NeoWs Live Test

Your challenge is to develop an interactive visualization and simulation tool that enables users to explore asteroid impact scenarios, predict consequences, and evaluate mitigation strategies using real NASA and USGS datasets.

Prerequisites:
- Python 3.9+
- `requests` and `python-dotenv` installed (for example: `pip3 install --user requests python-dotenv`)
- NASA API key stored in `.env` as `NASA_API_KEY=your_key`

Run the integration test (hits the live NASA NeoWs API):

```bash
python3 -m unittest tests.test_nasa_neo_live
```

If the environment variable is missing the live checks will be skipped automatically.
