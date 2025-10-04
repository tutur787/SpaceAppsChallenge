import os
import unittest

from src.api.nasa_neo import NeoWsError
from src.data.defaults import (
    DEFAULT_NEO_ID,
    DEFAULT_SBDB_ID,
    MATERIAL_PRESETS,
    SBDBError,
    fetch_neows_object,
    fetch_sbdb_payload,
    get_reference_defaults,
)
from tests.run_data_source_check import PARAMETERS as WHITEPAPER_PARAMS, get_nested


class TestDashboardDefaultsLive(unittest.TestCase):
    """Exercise the live data defaults that seed the Streamlit dashboard."""

    @classmethod
    def setUpClass(cls) -> None:
        if not os.getenv("NASA_API_KEY"):
            raise unittest.SkipTest("NASA_API_KEY not set; skip live dashboard default test")

    def test_reference_defaults_matches_live_data(self) -> None:
        try:
            defaults = get_reference_defaults()
        except (NeoWsError, SBDBError) as exc:
            self.skipTest(f"Live data unavailable: {exc}")

        neo_payload = fetch_neows_object(DEFAULT_NEO_ID)
        sbdb_payload = fetch_sbdb_payload(DEFAULT_SBDB_ID)

        print("\n=== NeoWs payload (top-level columns) ===")
        for key in sorted(neo_payload.keys()):
            print(f"  {key}")

        print("\n=== SBDB payload (top-level columns) ===")
        for key in sorted(sbdb_payload.keys()):
            print(f"  {key}")

        phys_par = sbdb_payload.get("phys_par")
        if isinstance(phys_par, dict):
            print("  phys_par ->")
            for key in sorted(phys_par.keys()):
                print(f"    {key}")

        print("\n=== White paper parameters & data sources ===")
        header = f"{ 'Parameter':35} {'NeoWs value':20} {'SBDB value':20} {'Dashboard':20}"
        print(header)
        print("-" * len(header))

        dashboard_field_map = {
            "Diameter max (m)": defaults.get("diameter_m"),
            "Bulk density (kg/m^3)": defaults.get("density"),
            "Relative velocity (km/s)": defaults.get("velocity_km_s"),
            "Spectral class": defaults.get("material"),
        }

        def _fmt(value: object) -> str:
            if value is None or value == "":
                return "—"
            if isinstance(value, float):
                return f"{value:,.3g}"
            if isinstance(value, (int,)):
                return f"{value:,.0f}"
            text = str(value)
            return text if len(text) <= 19 else text[:16] + "…"

        for spec in WHITEPAPER_PARAMS:
            neo_value = get_nested(neo_payload, spec.neo_path)
            sbdb_value = get_nested(sbdb_payload, spec.sbdb_path)
            dashboard_value = dashboard_field_map.get(spec.label, "not used")
            print(
                f"{spec.label:35} "
                f"{_fmt(neo_value):20} "
                f"{_fmt(sbdb_value):20} "
                f"{_fmt(dashboard_value):20}"
            )

        neo_name = defaults.get("neo_name") or "(unknown)"
        sbdb_name = defaults.get("sbdb_fullname") or DEFAULT_SBDB_ID
        diameter = defaults.get("diameter_m") or 0.0
        velocity = defaults.get("velocity_km_s") or 0.0
        density = defaults.get("density") or 0.0
        material = defaults.get("material") or "Stony (ordinary chondrite)"
        provenance = defaults.get("provenance", "synthetic fallback")

        print("\n=== Dashboard defaults (mirrors sidebar) ===")
        print(f"  NeoWs object      : {neo_name} (ID {DEFAULT_NEO_ID})")
        print(f"  SBDB object       : {sbdb_name}")
        print(f"  Material preset   : {material}")
        print(f"  Provenance        : {provenance}")
        print(f"  Diameter (m)      : {diameter:,.1f}")
        print(f"  Velocity (km/s)   : {velocity:,.2f}")
        print(f"  Bulk density (kg/m^3): {density:,.0f}")

        self.assertGreater(diameter, 0.0, "Expected live diameter > 0")
        self.assertGreater(velocity, 0.0, "Expected live velocity > 0")
        self.assertGreater(density, 0.0, "Expected live density > 0")
        self.assertIn(material, MATERIAL_PRESETS, "Material should be one of the known presets")
        self.assertNotEqual(provenance, "synthetic fallback", "Live defaults should not fall back to synthetic data")


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
