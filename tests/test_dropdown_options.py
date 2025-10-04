"""
Smoke tests for dropdown option validation.

Phase 3 acceptance checks: ensure dropdown options are non-empty and well-formed.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
from main import MATERIAL_PRESETS, CITY_PRESETS, CITY_RING_DENSITY


class TestMaterialPresetDropdown:
    """Tests for material preset dropdown options."""

    def test_material_presets_not_empty(self):
        """Material preset dropdown must have at least one option."""
        assert len(MATERIAL_PRESETS) > 0, "MATERIAL_PRESETS should not be empty"

    def test_material_presets_count(self):
        """Material preset dropdown should have exactly 4 options per roadmap."""
        assert len(MATERIAL_PRESETS) == 4, f"Expected 4 material presets, found {len(MATERIAL_PRESETS)}"

    def test_material_preset_keys_valid(self):
        """All material preset keys should be non-empty strings."""
        for key in MATERIAL_PRESETS.keys():
            assert isinstance(key, str), f"Material preset key {key!r} is not a string"
            assert len(key) > 0, "Material preset key should not be empty"

    def test_material_preset_structure(self):
        """Each material preset should have density and strength_mpa fields."""
        for name, preset in MATERIAL_PRESETS.items():
            assert "density" in preset, f"Material preset {name!r} missing 'density' field"
            assert "strength_mpa" in preset, f"Material preset {name!r} missing 'strength_mpa' field"
            assert isinstance(preset["density"], (int, float)), f"Material preset {name!r} density is not numeric"
            assert isinstance(preset["strength_mpa"], (int, float)), f"Material preset {name!r} strength_mpa is not numeric"
            assert preset["density"] > 0, f"Material preset {name!r} density should be positive"
            assert preset["strength_mpa"] > 0, f"Material preset {name!r} strength_mpa should be positive"


class TestCityPresetDropdown:
    """Tests for city preset dropdown options."""

    def test_city_presets_not_empty(self):
        """City preset dropdown must have at least one option."""
        assert len(CITY_PRESETS) > 0, "CITY_PRESETS should not be empty"

    def test_city_presets_count(self):
        """City preset dropdown should have exactly 7 options per roadmap."""
        assert len(CITY_PRESETS) == 7, f"Expected 7 city presets, found {len(CITY_PRESETS)}"

    def test_city_preset_keys_valid(self):
        """All city preset keys should be non-empty strings."""
        for key in CITY_PRESETS.keys():
            assert isinstance(key, str), f"City preset key {key!r} is not a string"
            assert len(key) > 0, "City preset key should not be empty"

    def test_city_preset_structure(self):
        """Each city preset should have (lat, lon) tuple or (None, None) for placeholder."""
        for name, coords in CITY_PRESETS.items():
            assert isinstance(coords, tuple), f"City preset {name!r} coordinates are not a tuple"
            assert len(coords) == 2, f"City preset {name!r} should have exactly 2 coordinates"
            lat, lon = coords
            if lat is not None and lon is not None:
                assert isinstance(lat, (int, float)), f"City preset {name!r} latitude is not numeric"
                assert isinstance(lon, (int, float)), f"City preset {name!r} longitude is not numeric"
                assert -90 <= lat <= 90, f"City preset {name!r} latitude {lat} out of range [-90, 90]"
                assert -180 <= lon <= 180, f"City preset {name!r} longitude {lon} out of range [-180, 180]"

    def test_city_preset_has_placeholder(self):
        """City presets should include a placeholder option."""
        placeholders = [name for name, coords in CITY_PRESETS.items() if coords == (None, None)]
        assert len(placeholders) > 0, "CITY_PRESETS should include at least one placeholder option"


class TestCityRingDensityMapping:
    """Tests for city ring density mappings."""

    def test_city_ring_density_not_empty(self):
        """City ring density mapping must have at least one entry."""
        assert len(CITY_RING_DENSITY) > 0, "CITY_RING_DENSITY should not be empty"

    def test_city_ring_density_count(self):
        """City ring density mapping should have exactly 6 cities per roadmap."""
        assert len(CITY_RING_DENSITY) == 6, f"Expected 6 city density mappings, found {len(CITY_RING_DENSITY)}"

    def test_city_ring_density_structure(self):
        """Each city ring density should have severe, moderate, and light fields."""
        for city, rings in CITY_RING_DENSITY.items():
            assert "severe" in rings, f"City {city!r} missing 'severe' ring density"
            assert "moderate" in rings, f"City {city!r} missing 'moderate' ring density"
            assert "light" in rings, f"City {city!r} missing 'light' ring density"
            for ring_name, density in rings.items():
                assert isinstance(density, (int, float)), f"City {city!r} ring {ring_name!r} density is not numeric"
                assert density > 0, f"City {city!r} ring {ring_name!r} density should be positive"

    def test_city_ring_density_aligns_with_presets(self):
        """City ring density entries should correspond to actual city presets (excluding placeholder)."""
        valid_cities = {name for name, coords in CITY_PRESETS.items() if coords != (None, None)}
        density_cities = set(CITY_RING_DENSITY.keys())
        assert density_cities.issubset(valid_cities), (
            f"CITY_RING_DENSITY contains cities not in CITY_PRESETS: "
            f"{density_cities - valid_cities}"
        )


class TestDropdownConsistency:
    """Cross-dropdown consistency checks."""

    def test_all_dropdowns_have_options(self):
        """All dropdown data structures should be non-empty."""
        assert len(MATERIAL_PRESETS) > 0, "MATERIAL_PRESETS is empty"
        assert len(CITY_PRESETS) > 0, "CITY_PRESETS is empty"
        assert len(CITY_RING_DENSITY) > 0, "CITY_RING_DENSITY is empty"

    def test_material_preset_names_unique(self):
        """Material preset names should be unique."""
        names = list(MATERIAL_PRESETS.keys())
        assert len(names) == len(set(names)), "MATERIAL_PRESETS has duplicate names"

    def test_city_preset_names_unique(self):
        """City preset names should be unique."""
        names = list(CITY_PRESETS.keys())
        assert len(names) == len(set(names)), "CITY_PRESETS has duplicate names"
