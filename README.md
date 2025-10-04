# SpaceAppsChallenge

## NeoWs Live Test

Prerequisites:
- Python 3.9+
- `requests` and `python-dotenv` installed (for example: `pip3 install --user requests python-dotenv`)
- NASA API key stored in `.env` as `NASA_API_KEY=your_key`

Run the integration test (hits the live NASA NeoWs API):

```bash
python3 -m unittest tests.test_nasa_neo_live
```

If the environment variable is missing the live checks will be skipped automatically.
