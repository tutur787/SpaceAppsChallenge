import json
from pathlib import Path

BASE = Path(__file__).parent.parent / "src" / "i18n"


def _flatten(d, prefix=""):
    out = {}
    for k, v in d.items():
        key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict):
            out.update(_flatten(v, key))
        else:
            out[key] = v
    return out


def test_translation_key_parity():
    en = json.load((BASE / "en.json").open("r", encoding="utf-8"))
    en_keys = set(_flatten(en).keys())

    for p in BASE.glob("*.json"):
        if p.name == "en.json":
            continue
        other = json.load(p.open("r", encoding="utf-8"))
        other_keys = set(_flatten(other).keys())
        missing = en_keys - other_keys
        assert not missing, f"Missing keys in {p.name}: {sorted(missing)[:10]}"
