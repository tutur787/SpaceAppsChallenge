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


def test_slider_telemetry_records_status_and_provenance():
    state = {}
    report = ensure_report(state)

    update_slider_validation(report, "diameter_m", "ok")
    update_slider_provenance(
        report,
        "diameter_m",
        ProvenanceTag.LIVE,
        value=120.0,
        detail="Hydrated from NeoWs",
        widget_key="diameter_m",
    )

    slider_entry = state["telemetry_report"]["sliders"]["diameter_m"]
    assert slider_entry["status"] == "ok"
    assert slider_entry["provenance"] == ProvenanceTag.LIVE.value
    assert slider_entry["value"] == 120.0
    assert slider_entry["widget_key"] == "diameter_m"


def test_dataset_provenance_and_formatting():
    state = {}
    report = ensure_report(state)

    record_dataset_provenance(
        report,
        "neo_feed",
        ProvenanceTag.LIVE,
        status="ok",
        detail="NeoWs feed hydrated.",
    )
    add_warning(report, "Fallback defaults in use.")

    rendered = format_report(state["telemetry_report"])
    assert "neo_feed" in rendered
    assert "Fallback defaults in use." in rendered


def test_synthetic_provenance_tags():
    state = {}
    report = ensure_report(state)

    # Test synthetic population density tracking
    record_dataset_provenance(
        report,
        "population_density",
        ProvenanceTag.SYNTHETIC_POP_DENSITY,
        status="warning",
        detail="Using DEFAULT_RING_DENSITY synthetic fallback.",
    )

    # Test synthetic casualty rate tracking
    record_dataset_provenance(
        report,
        "casualty_rates",
        ProvenanceTag.SYNTHETIC_CASUALTY,
        status="warning",
        detail="Using SYNTHETIC_CASUALTY_RATE fallback values.",
    )

    # Test city preset density tracking
    record_dataset_provenance(
        report,
        "population_density_city",
        ProvenanceTag.CITY_PRESET_DENSITY,
        status="ok",
        detail="Using city preset density for Tokyo, Japan.",
    )

    rendered = format_report(state["telemetry_report"])
    assert "synthetic_pop_density" in rendered
    assert "synthetic_casualty" in rendered
    assert "city_preset_density" in rendered
    assert "DEFAULT_RING_DENSITY" in rendered


def test_error_tracking_and_formatting():
    state = {}
    report = ensure_report(state)

    add_error(report, "Synthetic data used in live mode without approval.")
    add_warning(report, "NeoWs feed unavailable.")

    rendered = format_report(state["telemetry_report"])
    assert "Errors:" in rendered
    assert "Synthetic data used in live mode" in rendered
    assert "Warnings:" in rendered
    assert "NeoWs feed unavailable" in rendered


def test_live_mode_synthetic_data_validation():
    """Test that synthetic data usage in live mode is flagged."""
    state = {}
    report = ensure_report(state)

    # Simulate live mode with synthetic data
    record_dataset_provenance(
        report,
        "neo_feed",
        ProvenanceTag.LIVE,
        status="ok",
        detail="NeoWs daily feed for 2025-10-04.",
    )

    record_dataset_provenance(
        report,
        "population_density",
        ProvenanceTag.SYNTHETIC_POP_DENSITY,
        status="warning",
        detail="Using DEFAULT_RING_DENSITY synthetic fallback.",
    )

    add_warning(
        report,
        "Using synthetic population density in live mode — consider switching to a city preset or integrating real population data.",
    )

    datasets = state["telemetry_report"]["datasets"]
    assert datasets["neo_feed"]["provenance"] == "live"
    assert datasets["population_density"]["provenance"] == "synthetic_pop_density"
    assert datasets["population_density"]["status"] == "warning"

    warnings = state["telemetry_report"]["warnings"]
    assert any("synthetic population density in live mode" in w for w in warnings)


def test_dropdown_telemetry():
    """Test dropdown selection tracking."""
    state = {}
    report = ensure_report(state)

    # Track material preset dropdown
    record_dropdown_selection(
        report,
        "material_preset",
        "Stony (ordinary chondrite)",
        4,
        provenance=ProvenanceTag.MATERIAL_PRESET,
        detail="Selected material: Stony (ordinary chondrite)",
    )

    # Track city preset dropdown
    record_dropdown_selection(
        report,
        "city_preset",
        "Tokyo, Japan",
        7,
        provenance=ProvenanceTag.CITY_PRESET_DENSITY,
        detail="Selected location: Tokyo, Japan",
    )

    dropdowns = state["telemetry_report"]["dropdowns"]
    assert dropdowns["material_preset"]["selected"] == "Stony (ordinary chondrite)"
    assert dropdowns["material_preset"]["total_options"] == 4
    assert dropdowns["material_preset"]["provenance"] == "material_preset"

    assert dropdowns["city_preset"]["selected"] == "Tokyo, Japan"
    assert dropdowns["city_preset"]["total_options"] == 7
    assert dropdowns["city_preset"]["provenance"] == "city_preset_density"

    rendered = format_report(state["telemetry_report"])
    assert "Dropdowns:" in rendered
    assert "material_preset" in rendered
    assert "city_preset" in rendered
    assert "total_options=4" in rendered
    assert "total_options=7" in rendered
