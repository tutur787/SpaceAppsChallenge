from typing import Set

from src.config import SLIDER_SPECS

# Slider keys that must be represented in the canonical map to satisfy Phase 1 roadmap work.
EXPECTED_SLIDER_KEYS: Set[str] = {
    "diameter_m",
    "velocity_km_s",
    "angle_deg",
    "bulk_density",
    "strength_mpa",
    "pair_samples",
    "delta_v_mm_s",
    "lead_days",
    "inbound_bearing",
}


def test_slider_specs_cover_expected_keys():
    missing = EXPECTED_SLIDER_KEYS - set(SLIDER_SPECS)
    assert not missing, f"Missing slider specs for: {sorted(missing)}"


def test_slider_specs_have_valid_bounds():
    for key, spec in SLIDER_SPECS.items():
        assert spec.label, f"Spec '{key}' should define a label"
        assert spec.min_value <= spec.max_value, f"Spec '{key}' min > max"
        if spec.whitepaper_default is not None:
            assert spec.min_value <= spec.whitepaper_default <= spec.max_value, (
                f"Spec '{key}' default {spec.whitepaper_default} outside [{spec.min_value}, {spec.max_value}]"
            )
        assert spec.step > 0, f"Spec '{key}' step must be positive"
