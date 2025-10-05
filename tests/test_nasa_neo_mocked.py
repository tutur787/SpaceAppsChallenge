"""Mocked tests for NASA NeoWs client that don't require network access."""
import unittest
from unittest.mock import MagicMock, patch

import requests

from src.api.nasa_neo import (
    NeoWsClient,
    NeoWsError,
    estimate_mean_diameter_km,
    extract_orbital_elements,
)


class TestNeoWsClientMocked(unittest.TestCase):
    """Mocked tests for NeoWs client that work without network access."""

    @patch("src.api.nasa_neo.requests.Session")
    def test_get_feed_parses_response(self, mock_session_class: MagicMock) -> None:
        """Test that get_feed correctly parses a mocked API response."""
        mock_session = MagicMock()
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "near_earth_objects": {
                "2023-01-01": [
                    {
                        "name": "Test Asteroid",
                        "neo_reference_id": "123456",
                        "estimated_diameter": {
                            "meters": {
                                "estimated_diameter_min": 100.0,
                                "estimated_diameter_max": 200.0,
                            }
                        },
                        "is_potentially_hazardous_asteroid": True,
                        "close_approach_data": [
                            {
                                "relative_velocity": {"kilometers_per_second": "15.5"},
                                "miss_distance": {"kilometers": "1000000"},
                            }
                        ],
                    }
                ]
            }
        }
        mock_session.get.return_value = mock_response
        mock_session_class.return_value = mock_session

        with NeoWsClient(api_key="TEST_KEY") as client:
            payload = client.get_feed(start_date="2023-01-01", end_date="2023-01-02")

        self.assertIn("near_earth_objects", payload)
        self.assertIn("2023-01-01", payload["near_earth_objects"])
        self.assertEqual(len(payload["near_earth_objects"]["2023-01-01"]), 1)

    @patch("src.api.nasa_neo.requests.Session")
    def test_get_neo_parses_orbital_elements(self, mock_session_class: MagicMock) -> None:
        """Test that get_neo correctly parses orbital elements."""
        mock_session = MagicMock()
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "name": "(2023 ABC)",
            "neo_reference_id": "3542519",
            "orbital_data": {
                "semi_major_axis": "1.5",
                "eccentricity": "0.2",
                "inclination": "5.5",
                "ascending_node_longitude": "120.0",
                "perihelion_argument": "90.0",
                "mean_anomaly": "45.0",
                "epoch_osculation": "2459000.5",
            },
        }
        mock_session.get.return_value = mock_response
        mock_session_class.return_value = mock_session

        with NeoWsClient(api_key="TEST_KEY") as client:
            payload = client.get_neo("3542519")

        self.assertEqual(payload["name"], "(2023 ABC)")
        elements = extract_orbital_elements(payload)
        self.assertEqual(elements["semi_major_axis_au"], 1.5)
        self.assertEqual(elements["eccentricity"], 0.2)
        self.assertEqual(elements["inclination_deg"], 5.5)

    @patch("src.api.nasa_neo.requests.Session")
    def test_network_error_raises_neows_error(self, mock_session_class: MagicMock) -> None:
        """Test that network errors are wrapped in NeoWsError."""
        mock_session = MagicMock()
        mock_session.get.side_effect = requests.RequestException("Network error")
        mock_session_class.return_value = mock_session

        with self.assertRaises(NeoWsError):
            with NeoWsClient(api_key="TEST_KEY") as client:
                client.get_feed(start_date="2023-01-01")

    def test_extract_orbital_elements_handles_missing_data(self) -> None:
        """Test that orbital element extraction handles missing data gracefully."""
        payload = {"name": "Test", "orbital_data": {}}
        elements = extract_orbital_elements(payload)
        self.assertIsNone(elements["semi_major_axis_au"])
        self.assertIsNone(elements["eccentricity"])

    def test_estimate_diameter_computes_average(self) -> None:
        """Test diameter estimation computes the mean correctly."""
        payload = {
            "estimated_diameter": {
                "kilometers": {
                    "estimated_diameter_min": 0.1,
                    "estimated_diameter_max": 0.3,
                }
            }
        }
        diameter = estimate_mean_diameter_km(payload)
        self.assertAlmostEqual(diameter, 0.2)

    def test_estimate_diameter_handles_missing_data(self) -> None:
        """Test diameter estimation returns None for missing data."""
        payload = {"estimated_diameter": {}}
        diameter = estimate_mean_diameter_km(payload)
        self.assertIsNone(diameter)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
