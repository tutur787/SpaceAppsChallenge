import os
import unittest

from src.api.nasa_neo import NeoWsClient, NeoWsError, extract_orbital_elements


class TestNeoWsClientLive(unittest.TestCase):
    """Live tests that reach out to NASA NeoWs. Requires NASA_API_KEY."""

    @classmethod
    def setUpClass(cls) -> None:
        if not os.getenv("NASA_API_KEY"):
            raise unittest.SkipTest("NASA_API_KEY not set; skip live NeoWs tests")

    def test_get_feed_returns_objects(self) -> None:
        with NeoWsClient() as client:
            payload = client.get_feed(start_date="2023-01-01", end_date="2023-01-02")
        self.assertIn("near_earth_objects", payload)
        self.assertTrue(payload["near_earth_objects"], "Expected NEO data in feed response")

    def test_get_neo_returns_orbital_elements(self) -> None:
        with NeoWsClient() as client:
            payload = client.get_neo("3542519")
        self.assertIn("name", payload)
        self.assertTrue(payload["name"])
        elements = extract_orbital_elements(payload)
        self.assertTrue(
            any(value is not None for value in elements.values()),
            "Expected at least one orbital element to be parsed",
        )


class TestNeoWsKeyHandling(unittest.TestCase):
    def test_missing_key_raises_error(self) -> None:
        original = os.environ.get("NASA_API_KEY")
        try:
            if "NASA_API_KEY" in os.environ:
                del os.environ["NASA_API_KEY"]
            with self.assertRaises(NeoWsError):
                NeoWsClient(api_key=None)
        finally:
            if original is not None:
                os.environ["NASA_API_KEY"] = original


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
