"""Regression tests for headless telemetry output.

These tests verify that:
1. Telemetry reports format correctly
2. Exit codes are non-zero when errors exist
3. Reports include all expected sections (datasets, sliders, dropdowns, warnings, errors)
"""

from src.telemetry import (
    ProvenanceTag,
    add_error,
    add_warning,
    ensure_report,
    format_report,
    record_dataset_provenance,
    record_dropdown_selection,
    update_slider_provenance,
    update_slider_validation,
)


def test_format_report_empty():
    """Empty report should produce minimal output."""
    report = {"sliders": {}, "datasets": {}, "warnings": []}
    output = format_report(report)
    assert "Telemetry report" in output
    assert "Datasets:" not in output
    assert "Sliders:" not in output
    assert "Warnings:" not in output


def test_format_report_with_datasets():
    """Report should include dataset section when datasets exist."""
    report = {"sliders": {}, "datasets": {}, "warnings": []}
    record_dataset_provenance(
        report, "neo_data", ProvenanceTag.LIVE, status="ok", detail="Fetched from NeoWs"
    )

    output = format_report(report)
    assert "Datasets:" in output
    assert "neo_data" in output
    assert "status=ok" in output
    assert "provenance=live" in output
    assert "Fetched from NeoWs" in output


def test_format_report_with_sliders():
    """Report should include slider section when sliders exist."""
    report = {"sliders": {}, "datasets": {}, "warnings": []}
    update_slider_provenance(
        report,
        "diameter",
        ProvenanceTag.LIVE,
        value=100.0,
        detail="From selected NEO",
    )
    update_slider_validation(report, "diameter", "valid")

    output = format_report(report)
    assert "Sliders:" in output
    assert "diameter" in output
    assert "status=valid" in output
    assert "provenance=live" in output
    assert "value=100.0" in output
    assert "From selected NEO" in output


def test_format_report_with_dropdowns():
    """Report should include dropdown section when dropdowns exist."""
    report = {"sliders": {}, "datasets": {}, "warnings": []}
    record_dropdown_selection(
        report,
        "material_preset",
        "Stony (ordinary chondrite)",
        total_options=4,
        provenance=ProvenanceTag.MATERIAL_PRESET,
    )

    output = format_report(report)
    assert "Dropdowns:" in output
    assert "material_preset" in output
    assert "selected='Stony (ordinary chondrite)'" in output
    assert "total_options=4" in output
    assert "provenance=material_preset" in output


def test_format_report_with_warnings():
    """Report should include warnings section."""
    report = {"sliders": {}, "datasets": {}, "warnings": []}
    add_warning(report, "Synthetic data used for population density")

    output = format_report(report)
    assert "Warnings:" in output
    assert "Synthetic data used for population density" in output


def test_format_report_with_errors():
    """Report should include errors section."""
    report = {"sliders": {}, "datasets": {}, "warnings": [], "errors": []}
    add_error(report, "Failed to fetch live NEO data")

    output = format_report(report)
    assert "Errors:" in output
    assert "Failed to fetch live NEO data" in output


def test_format_report_comprehensive():
    """Report should include all sections when populated."""
    report = {"sliders": {}, "datasets": {}, "warnings": [], "errors": []}

    # Add dataset
    record_dataset_provenance(
        report, "neo_data", ProvenanceTag.LIVE, status="ok"
    )

    # Add slider
    update_slider_provenance(
        report, "diameter", ProvenanceTag.LIVE, value=100.0
    )
    update_slider_validation(report, "diameter", "valid")

    # Add dropdown
    record_dropdown_selection(
        report, "city_preset", "New York City, USA", total_options=7
    )

    # Add warning
    add_warning(report, "Using approximate crater model")

    # Add error
    add_error(report, "Validation failed for angle slider")

    output = format_report(report)

    # Verify all sections present
    assert "Datasets:" in output
    assert "neo_data" in output

    assert "Sliders:" in output
    assert "diameter" in output

    assert "Dropdowns:" in output
    assert "city_preset" in output

    assert "Warnings:" in output
    assert "Using approximate crater model" in output

    assert "Errors:" in output
    assert "Validation failed for angle slider" in output


def test_healthy_run_has_no_errors():
    """A healthy telemetry run should have no errors."""
    report = ensure_report({})

    # Populate with valid data
    record_dataset_provenance(report, "neo_data", ProvenanceTag.LIVE, status="ok")
    update_slider_provenance(report, "diameter", ProvenanceTag.LIVE, value=100.0)
    update_slider_validation(report, "diameter", "valid")

    # Verify no errors
    assert report.get("errors", []) == []


def test_error_exit_code_detection():
    """Test that errors are detectable for exit code purposes."""
    report = ensure_report({})

    # Initially no errors
    assert len(report.get("errors", [])) == 0

    # Add an error
    add_error(report, "Critical validation failure")

    # Now errors exist
    errors = report.get("errors", [])
    assert len(errors) == 1
    assert errors[0] == "Critical validation failure"


def test_telemetry_sections_ordered():
    """Verify telemetry report sections appear in expected order."""
    report = {"sliders": {}, "datasets": {}, "warnings": [], "errors": []}

    record_dataset_provenance(report, "neo", ProvenanceTag.LIVE)
    update_slider_provenance(report, "vel", ProvenanceTag.LIVE, value=20.0)
    record_dropdown_selection(report, "mat", "Iron-nickel", total_options=4)
    add_warning(report, "Test warning")
    add_error(report, "Test error")

    output = format_report(report)
    lines = output.split("\n")

    # Find section indices
    datasets_idx = next(i for i, line in enumerate(lines) if "Datasets:" in line)
    sliders_idx = next(i for i, line in enumerate(lines) if "Sliders:" in line)
    dropdowns_idx = next(i for i, line in enumerate(lines) if "Dropdowns:" in line)
    warnings_idx = next(i for i, line in enumerate(lines) if "Warnings:" in line)
    errors_idx = next(i for i, line in enumerate(lines) if "Errors:" in line)

    # Verify ordering
    assert datasets_idx < sliders_idx
    assert sliders_idx < dropdowns_idx
    assert dropdowns_idx < warnings_idx
    assert warnings_idx < errors_idx
