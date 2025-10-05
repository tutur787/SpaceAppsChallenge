#!/usr/bin/env bash
set -e
# Use the port Azure provides, default to 8000 if missing
export PORT="${PORT:-8000}"
python -m streamlit run main.py --server.port="$PORT" --server.address=0.0.0.0 --browser.gatherUsageStats=false
