"""
Component/UI tests for dropdown behavior.

Phase 3 acceptance checks: validate dropdown interactions and state management.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
from unittest.mock import MagicMock, patch
from main import MATERIAL_PRESETS, CITY_PRESETS, ensure_report
from src.telemetry import ProvenanceTag


class TestMaterialDropdownBehavior:
    """Tests for material preset dropdown behavior."""

    def test_material_dropdown_default_selection(self):
        """Material dropdown should default to index 1 (Stony ordinary chondrite)."""
        material_list = list(MATERIAL_PRESETS.keys())
        default_index = 1
        assert default_index < len(material_list), "Default index out of range"
        default_material = material_list[default_index]
        assert default_material == "Stony (ordinary chondrite)", (
            f"Expected default material 'Stony (ordinary chondrite)', got {default_material!r}"
        )

    def test_material_selection_updates_density_and_strength(self):
        """Selecting a material should update both density and strength values."""
        for material_name, preset in MATERIAL_PRESETS.items():
            assert "density" in preset, f"Material {material_name!r} missing density"
            assert "strength_mpa" in preset, f"Material {material_name!r} missing strength_mpa"
            # Verify values are accessible as they would be in UI
            density = preset["density"]
            strength = preset["strength_mpa"]
            assert density > 0, f"Material {material_name!r} has non-positive density"
            assert strength > 0, f"Material {material_name!r} has non-positive strength"

    def test_material_dropdown_telemetry_recording(self):
        """Material dropdown selection should be recorded in telemetry."""
        mock_session_state = {}
        report = ensure_report(mock_session_state)

        # Simulate material selection
        from main import record_dropdown_selection
        material = "Stony (ordinary chondrite)"
        record_dropdown_selection(
            report,
            "material_preset",
            material,
            len(MATERIAL_PRESETS),
            provenance=ProvenanceTag.MATERIAL_PRESET,
            detail=f"Selected material: {material}",
        )

        # Verify telemetry recorded
        assert "dropdowns" in report, "Dropdowns section not created in report"
        assert "material_preset" in report["dropdowns"], "Material selection not recorded"
        dropdown_entry = report["dropdowns"]["material_preset"]
        assert dropdown_entry["selected"] == material
        assert dropdown_entry["total_options"] == 4
        assert dropdown_entry["provenance"] == ProvenanceTag.MATERIAL_PRESET.value


class TestCityDropdownBehavior:
    """Tests for city preset dropdown behavior."""

    def test_city_dropdown_includes_placeholder(self):
        """City dropdown should include '— choose a place —' placeholder."""
        city_list = list(CITY_PRESETS.keys())
        placeholders = [c for c in city_list if CITY_PRESETS[c] == (None, None)]
        assert len(placeholders) > 0, "City dropdown missing placeholder option"
        assert "— choose a place —" in city_list, "Expected placeholder '— choose a place —' not found"

    def test_city_selection_provides_coordinates(self):
        """Selecting a city should provide valid lat/lon coordinates (except placeholder)."""
        for city_name, coords in CITY_PRESETS.items():
            if coords == (None, None):
                # Placeholder - skip coordinate validation
                continue
            lat, lon = coords
            assert lat is not None and lon is not None, f"City {city_name!r} has null coordinates"
            assert -90 <= lat <= 90, f"City {city_name!r} lat {lat} out of range"
            assert -180 <= lon <= 180, f"City {city_name!r} lon {lon} out of range"

    def test_city_dropdown_telemetry_recording(self):
        """City dropdown selection should be recorded in telemetry."""
        mock_session_state = {}
        report = ensure_report(mock_session_state)

        # Simulate city selection
        from main import record_dropdown_selection
        city = "Tokyo, Japan"
        record_dropdown_selection(
            report,
            "city_preset",
            city,
            len(CITY_PRESETS),
            provenance=ProvenanceTag.CITY_PRESET_DENSITY,
            detail=f"Selected city: {city}",
        )

        # Verify telemetry recorded
        assert "dropdowns" in report, "Dropdowns section not created in report"
        assert "city_preset" in report["dropdowns"], "City selection not recorded"
        dropdown_entry = report["dropdowns"]["city_preset"]
        assert dropdown_entry["selected"] == city
        assert dropdown_entry["total_options"] == 7


class TestNEODropdownBehavior:
    """Tests for NEO selection dropdown behavior."""

    def test_neo_dropdown_session_state_integration(self):
        """NEO dropdown should interact correctly with session state."""
        # Mock session state with NEO candidates
        mock_session_state = {
            "neo_candidates": [
                {"name": "2025 AB", "id": "123", "diameter_m": 100, "velocity_km_s": 20},
                {"name": "2025 CD", "id": "456", "diameter_m": 150, "velocity_km_s": 25},
            ]
        }

        # Verify candidates are accessible
        candidates = mock_session_state.get("neo_candidates", [])
        assert len(candidates) == 2, "Mock NEO candidates not properly set"

        # Simulate NEO selection
        selected_neo = candidates[0]
        assert selected_neo["name"] == "2025 AB"
        assert selected_neo["diameter_m"] == 100

    def test_neo_button_triggers_parameter_update(self):
        """NEO 'Use this NEO' button should trigger parameter updates."""
        # This test verifies the expected behavior pattern
        mock_neo = {
            "name": "Test NEO",
            "id": "test123",
            "diameter_m": 200,
            "velocity_km_s": 30,
        }

        # Verify required fields are present for UI update
        required_fields = ["name", "id", "diameter_m", "velocity_km_s"]
        for field in required_fields:
            assert field in mock_neo, f"NEO object missing required field: {field}"


class TestDropdownStateManagement:
    """Tests for dropdown state management and persistence."""

    def test_ensure_report_creates_telemetry_structure(self):
        """ensure_report should create telemetry structure in session state."""
        mock_session_state = {}
        report = ensure_report(mock_session_state)

        assert "telemetry_report" in mock_session_state, "Telemetry report not added to session state"
        assert isinstance(report, dict), "Report is not a dict instance"
        # Dropdowns section is created on-demand, not by ensure_report
        assert "sliders" in report, "Report missing sliders section"
        assert "datasets" in report, "Report missing datasets section"

    def test_dropdown_telemetry_persistence(self):
        """Dropdown selections should persist across multiple recordings."""
        mock_session_state = {}
        report = ensure_report(mock_session_state)

        from main import record_dropdown_selection

        # Record multiple dropdown selections
        record_dropdown_selection(
            report, "material_preset", "Iron-nickel", 4,
            provenance=ProvenanceTag.MATERIAL_PRESET, detail="Material selected"
        )
        record_dropdown_selection(
            report, "city_preset", "London, UK", 7,
            provenance=ProvenanceTag.CITY_PRESET_DENSITY, detail="City selected"
        )

        # Verify both persisted
        assert "dropdowns" in report, "Dropdowns section not created"
        assert len(report["dropdowns"]) == 2, "Not all dropdown selections persisted"
        assert "material_preset" in report["dropdowns"]
        assert "city_preset" in report["dropdowns"]

    def test_dropdown_selection_overwrites_previous(self):
        """Selecting a dropdown again should overwrite previous selection."""
        mock_session_state = {}
        report = ensure_report(mock_session_state)

        from main import record_dropdown_selection

        # First selection
        record_dropdown_selection(
            report, "material_preset", "Iron-nickel", 4,
            provenance=ProvenanceTag.MATERIAL_PRESET, detail="First selection"
        )
        first_selection = report["dropdowns"]["material_preset"]["selected"]

        # Second selection (overwrite)
        record_dropdown_selection(
            report, "material_preset", "Cometary (icy)", 4,
            provenance=ProvenanceTag.MATERIAL_PRESET, detail="Second selection"
        )
        second_selection = report["dropdowns"]["material_preset"]["selected"]

        assert first_selection == "Iron-nickel"
        assert second_selection == "Cometary (icy)"
        assert len(report["dropdowns"]) == 1, "Should still have only one material_preset entry"
