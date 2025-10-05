"""Lightweight helpers for accumulating telemetry and provenance metadata."""

from __future__ import annotations

from enum import Enum
from typing import Dict, List, MutableMapping, Optional


class ProvenanceTag(str, Enum):
    """Canonical provenance labels used across the app."""

    LIVE = "live"
    CACHED = "cached"
    WHITEPAPER = "whitepaper"
    MATERIAL_PRESET = "material_preset"
    SESSION_STATE = "session_state"
    FALLBACK = "fallback_value"
    SPEC_FLOOR = "spec_floor"
    SYNTHETIC = "synthetic"
    SYNTHETIC_POP_DENSITY = "synthetic_pop_density"
    SYNTHETIC_CASUALTY = "synthetic_casualty"
    CITY_PRESET_DENSITY = "city_preset_density"
    SBDB = "sbdb"
    UNKNOWN = "unknown"


def ensure_report(container: MutableMapping[str, object]) -> Dict[str, object]:
    """Return the telemetry report stored in ``container`` (creating it if missing)."""

    report = container.setdefault(
        "telemetry_report",
        {"sliders": {}, "datasets": {}, "warnings": []},
    )
    return report  # type: ignore[return-value]


def update_slider_validation(
    report: Dict[str, object],
    key: str,
    status: str,
    *,
    message: Optional[str] = None,
) -> None:
    sliders: Dict[str, Dict[str, object]] = report.setdefault("sliders", {})  # type: ignore[assignment]
    entry = sliders.setdefault(key, {})
    entry.update({"status": status, "message": message})


def update_slider_provenance(
    report: Dict[str, object],
    key: str,
    provenance: ProvenanceTag,
    *,
    value: Optional[float] = None,
    detail: Optional[str] = None,
    widget_key: Optional[str] = None,
) -> None:
    sliders: Dict[str, Dict[str, object]] = report.setdefault("sliders", {})  # type: ignore[assignment]
    entry = sliders.setdefault(key, {})
    entry.update(
        {
            "provenance": provenance.value,
            "provenance_detail": detail,
            "value": value,
            "widget_key": widget_key or entry.get("widget_key"),
        }
    )


def record_dataset_provenance(
    report: Dict[str, object],
    dataset_key: str,
    provenance: ProvenanceTag,
    *,
    status: str = "ok",
    detail: Optional[str] = None,
) -> None:
    datasets: Dict[str, Dict[str, object]] = report.setdefault("datasets", {})  # type: ignore[assignment]
    entry = datasets.setdefault(dataset_key, {})
    entry.update({
        "provenance": provenance.value,
        "status": status,
    })
    if detail:
        entry["detail"] = detail


def add_warning(report: Dict[str, object], message: str) -> None:
    warnings: List[str] = report.setdefault("warnings", [])  # type: ignore[assignment]
    if message not in warnings:
        warnings.append(message)


def add_error(report: Dict[str, object], message: str) -> None:
    """Add an error message to the report (used for validation failures)."""
    errors: List[str] = report.setdefault("errors", [])  # type: ignore[assignment]
    if message not in errors:
        errors.append(message)


def record_dropdown_selection(
    report: Dict[str, object],
    dropdown_key: str,
    selected_value: str,
    total_options: int,
    *,
    provenance: Optional[ProvenanceTag] = None,
    detail: Optional[str] = None,
) -> None:
    """Record a dropdown selection for telemetry tracking."""
    dropdowns: Dict[str, Dict[str, object]] = report.setdefault("dropdowns", {})  # type: ignore[assignment]
    entry = dropdowns.setdefault(dropdown_key, {})
    entry.update({
        "selected": selected_value,
        "total_options": total_options,
    })
    if provenance:
        entry["provenance"] = provenance.value
    if detail:
        entry["detail"] = detail


def format_report(report: Dict[str, object]) -> str:
    """Render a human-readable summary suitable for headless logging."""

    lines: List[str] = ["Telemetry report"]

    datasets: Dict[str, Dict[str, object]] = report.get("datasets", {})  # type: ignore[assignment]
    if datasets:
        lines.append("Datasets:")
        for name in sorted(datasets):
            data = datasets[name]
            detail = f" — {data['detail']}" if data.get("detail") else ""
            lines.append(
                f"  - {name}: status={data.get('status', 'unknown')} provenance={data.get('provenance', 'unknown')}{detail}"
            )

    sliders: Dict[str, Dict[str, object]] = report.get("sliders", {})  # type: ignore[assignment]
    if sliders:
        lines.append("Sliders:")
        for key in sorted(sliders):
            data = sliders[key]
            detail = f" ({data['provenance_detail']})" if data.get("provenance_detail") else ""
            value = data.get("value")
            value_str = f", value={value}" if value is not None else ""
            status = data.get("status", "n/a")
            lines.append(
                f"  - {key}: status={status}, provenance={data.get('provenance', 'unknown')}{value_str}{detail}"
            )
            if data.get("message"):
                lines.append(f"    message: {data['message']}")

    dropdowns: Dict[str, Dict[str, object]] = report.get("dropdowns", {})  # type: ignore[assignment]
    if dropdowns:
        lines.append("Dropdowns:")
        for key in sorted(dropdowns):
            data = dropdowns[key]
            selected = data.get("selected", "unknown")
            total = data.get("total_options", 0)
            provenance = data.get("provenance", "unknown")
            lines.append(f"  - {key}: selected='{selected}', total_options={total}, provenance={provenance}")
            detail = data.get("detail")
            if detail:
                lines.append(f"    detail: {detail}")

    warnings: List[str] = report.get("warnings", [])  # type: ignore[assignment]
    if warnings:
        lines.append("Warnings:")
        for warning in warnings:
            lines.append(f"  - {warning}")

    errors: List[str] = report.get("errors", [])  # type: ignore[assignment]
    if errors:
        lines.append("Errors:")
        for error in errors:
            lines.append(f"  - {error}")

    return "\n".join(lines)

