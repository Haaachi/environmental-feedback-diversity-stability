"""YAML configuration loading with single-level extends inheritance."""

from __future__ import annotations

import copy
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[3]


def deep_update(base: dict, override: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_update(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load_config(path: str | Path) -> dict:
    path = Path(path)
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    with open(path) as f:
        cfg = yaml.safe_load(f)

    if "extends" in cfg:
        parent_path = path.parent / cfg.pop("extends")
        parent = load_config(parent_path)
        cfg = deep_update(parent, cfg)

    cfg["_config_path"] = str(path)
    return cfg


def require(cfg: dict, *keys: str) -> None:
    for key in keys:
        parts = key.split(".")
        node = cfg
        for p in parts:
            if not isinstance(node, dict) or p not in node:
                raise KeyError(f"Missing required config key: {key}")
            node = node[p]
