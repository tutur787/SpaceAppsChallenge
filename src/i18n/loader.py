"""Minimal i18n loader for Streamlit app.

Features:
- Load JSON translation files from src/i18n/*.json
- Provide t(key, **kwargs) for template interpolation using str.format
- set/get language and a simple cache
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, Optional

_ROOT = Path(__file__).parent

_CACHE: Dict[str, Dict[str, str]] = {}
_LANG_LABEL_CACHE: Dict[str, str] = {}

AVAILABLE_LANGS = []
_DEFAULT_LANG = "en"
_current_lang = _DEFAULT_LANG


def _load_lang(lang: str) -> Dict[str, str]:
    if lang in _CACHE:
        return _CACHE[lang]
    path = _ROOT / f"{lang}.json"
    if not path.exists():
        raise FileNotFoundError(f"translation file not found: {path}")
    with path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)
    # flatten nested dicts to dot.keys
    def _flatten(d: Dict[str, Any], prefix: str = "") -> Dict[str, str]:
        out: Dict[str, str] = {}
        for k, v in d.items():
            key = f"{prefix}.{k}" if prefix else k
            if isinstance(v, dict):
                out.update(_flatten(v, key))
            else:
                out[key] = str(v)
        return out

    flat = _flatten(data)
    _CACHE[lang] = flat
    return flat


def _discover_available() -> None:
    global AVAILABLE_LANGS
    files = [p.stem for p in _ROOT.glob("*.json")]
    AVAILABLE_LANGS = sorted(files)


def set_lang(lang: str) -> None:
    global _current_lang
    if not AVAILABLE_LANGS:
        _discover_available()
    if lang not in AVAILABLE_LANGS:
        raise ValueError(f"language not available: {lang}")
    _current_lang = lang


def get_lang() -> str:
    if not AVAILABLE_LANGS:
        _discover_available()
    return _current_lang


def get_language_label(lang: str, *, fallback: Optional[str] = None) -> str:
    """Return a human-readable label for the given language code.

    Falls back to `meta.language_native`, then `meta.language_english`, then the
    provided fallback or the language code itself.
    """
    if not AVAILABLE_LANGS:
        _discover_available()
    if lang in _LANG_LABEL_CACHE:
        return _LANG_LABEL_CACHE[lang]

    try:
        data = _load_lang(lang)
    except FileNotFoundError:
        label = fallback or lang
    else:
        label = (
            data.get("meta.language_native")
            or data.get("meta.language_english")
            or fallback
            or lang
        )

    _LANG_LABEL_CACHE[lang] = label
    return label


def t(key: str, default: Optional[str] = None, **kwargs: Any) -> str:
    """Translate a dotted key and interpolate with kwargs.

    Falls back to default then to the English translation, finally to key.
    """
    if not AVAILABLE_LANGS:
        _discover_available()
    lang = get_lang()
    try:
        data = _load_lang(lang)
    except FileNotFoundError:
        data = {}

    val = data.get(key)
    if val is None and lang != _DEFAULT_LANG:
        # fallback to default lang
        try:
            val = _load_lang(_DEFAULT_LANG).get(key)
        except FileNotFoundError:
            val = None
    if val is None:
        val = default if default is not None else key
    try:
        return val.format(**kwargs)
    except Exception:
        return val


# Populate available languages at import time
_discover_available()
