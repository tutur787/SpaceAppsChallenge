# conftest.py - pytest configuration that runs before all tests
import os
from pathlib import Path

def pytest_configure(config):
    """Load .env file before running tests"""
    try:
        from dotenv import load_dotenv
        env_path = Path(__file__).parent / ".env"
        if env_path.exists():
            load_dotenv(env_path)
            print(f"✅ Loaded environment variables from {env_path}")
    except ImportError:
        print("⚠️  python-dotenv not installed, .env file not loaded")
