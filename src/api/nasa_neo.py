from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any, Dict, Optional

import requests

try:
    from dotenv import load_dotenv
except ImportError:  # optional dependency
    load_dotenv = None
else:
    load_dotenv()


DEFAULT_BASE_URL = "https://api.nasa.gov/neo/rest/v1"


class NeoWsError(RuntimeError):
    """Raised when the NeoWs API returns an HTTP or parsing error."""


@dataclass
class NeoWsClient:
    api_key: Optional[str] = None
    base_url: str = DEFAULT_BASE_URL
    timeout: int = 30
    session: Optional[requests.Session] = None

    def __post_init__(self) -> None:
        if self.api_key is None:
            self.api_key = (
                os.getenv("NASA_API_KEY")
                or os.getenv("NEO_API_KEY")
                or os.getenv("NEOWS_API_KEY")
            )
        if not self.api_key:
            raise NeoWsError(
                "NASA API key missing. Set NASA_API_KEY in .env or pass api_key explicitly."
            )
        self._session: requests.Session = self.session or requests.Session()

    def close(self) -> None:
        self._session.close()

    def __enter__(self) -> NeoWsClient:
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        self.close()
        return False

    def get_feed(
        self,
        start_date: str,
        end_date: Optional[str] = None,
        detailed: bool = False,
    ) -> Dict[str, Any]:
        """Fetch close-approach data for a date range (YYYY-MM-DD)."""
        params: Dict[str, Any] = {"start_date": start_date}
        if end_date:
            params["end_date"] = end_date
        if detailed:
            params["detailed"] = "true"
        return self._request("feed", params)

    def get_neo(self, asteroid_id: str) -> Dict[str, Any]:
        """Fetch detailed orbital and physical properties for a specific NEO."""
        return self._request(f"neo/{asteroid_id}")

    def browse(self, page: int = 0, size: int = 20) -> Dict[str, Any]:
        """Iterate through the catalog of NEOs (useful for seeding local caches)."""
        params = {"page": page, "size": size}
        return self._request("neo/browse", params)

    def _request(self, path: str, params: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        query = dict(params or {})
        query.setdefault("api_key", self.api_key)
        url = self._build_url(path)
        try:
            response = self._session.get(url, params=query, timeout=self.timeout)
            response.raise_for_status()
        except requests.RequestException as exc:
            raise NeoWsError(f"NeoWs request failed: {exc}") from exc
        try:
            return response.json()
        except ValueError as exc:
            raise NeoWsError("NeoWs response was not valid JSON") from exc

    def _build_url(self, path: str) -> str:
        return f"{self.base_url.rstrip('/')}/{path.lstrip('/')}"


def extract_orbital_elements(payload: Dict[str, Any]) -> Dict[str, Optional[float]]:
    """Normalise orbital elements into floats (AU, degrees, epoch in JD)."""
    orbital = payload.get("orbital_data") or {}
    return {
        "semi_major_axis_au": _to_float(orbital.get("semi_major_axis")),
        "eccentricity": _to_float(orbital.get("eccentricity")),
        "inclination_deg": _to_float(orbital.get("inclination")),
        "ascending_node_longitude_deg": _to_float(orbital.get("ascending_node_longitude")),
        "argument_of_periapsis_deg": _to_float(orbital.get("perihelion_argument")),
        "mean_anomaly_deg": _to_float(orbital.get("mean_anomaly")),
        "epoch_jd": _to_float(orbital.get("epoch_osculation")),
    }


def estimate_mean_diameter_km(payload: Dict[str, Any]) -> Optional[float]:
    """Return the average estimated diameter in kilometres if present."""
    diameter = (payload.get("estimated_diameter") or {}).get("kilometers") or {}
    minimum = _to_float(diameter.get("estimated_diameter_min"))
    maximum = _to_float(diameter.get("estimated_diameter_max"))
    if minimum is None or maximum is None:
        return None
    return (minimum + maximum) / 2.0


def _to_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
