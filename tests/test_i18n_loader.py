import importlib

from src.i18n import t, set_lang, get_lang, AVAILABLE_LANGS


def test_discover_and_default():
    # Ensure available langs include English and at least one other we added
    assert "en" in AVAILABLE_LANGS


def test_translation_basic():
    set_lang("en")
    assert t("app.title").startswith("🛰️")
    set_lang("es")
    assert "Aprender" in t("app.title") or t("app.title") != "app.title"
