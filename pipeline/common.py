"""
common.py — shared helpers for the coral multiomics pipeline.

Everything that depends on the dataset (paths, sample-name decoding, thresholds)
is read from config.yaml so the analysis scripts stay generic.
"""
import os
import json
import yaml
import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
def load_config(path=None):
    """Load config.yaml. Falls back to config/config.yaml relative to repo root."""
    if path is None:
        here = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(here, "..", "config", "config.yaml")
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    # resolve paths relative to the repo root (parent of pipeline/)
    cfg["_root"] = os.path.abspath(os.path.join(os.path.dirname(path), ".."))
    return cfg


def rooted(cfg, relpath):
    """Resolve a config-relative path against the repo root."""
    if os.path.isabs(relpath):
        return relpath
    return os.path.join(cfg["_root"], relpath)


def ensure_dirs(cfg):
    for key in ("dir", "tables_dir", "figures_dir"):
        os.makedirs(rooted(cfg, cfg["output"][key]), exist_ok=True)
    os.makedirs(rooted(cfg, "results/_work"), exist_ok=True)


def work_path(cfg, name):
    """Path for intermediate artifacts shared between steps."""
    return os.path.join(rooted(cfg, "results/_work"), name)


def table_path(cfg, name):
    return os.path.join(rooted(cfg, cfg["output"]["tables_dir"]), name)


def figure_path(cfg, name):
    return os.path.join(rooted(cfg, cfg["output"]["figures_dir"]), name)


# --------------------------------------------------------------------------
def parse_sample(col, cfg):
    """Decode one sample column header into (timepoint_label, treatment).

    Driven entirely by the `sample_naming` block in config.yaml, so adapting to
    a new naming scheme only requires editing the config.
    """
    sn = cfg["sample_naming"]
    parts = col.split(sn["delimiter"])
    time_tok = parts[sn["time_token_index"]]
    cond_tok = parts[sn["cond_token_index"]]
    if cond_tok.startswith(sn["heat_prefix"]):
        treatment = "Heat"
    elif cond_tok.startswith(sn["control_prefix"]):
        treatment = "Control"
    else:
        raise ValueError(f"Cannot classify treatment for sample '{col}' "
                         f"(token '{cond_tok}')")
    time_label = sn["time_map"].get(time_tok)
    if time_label is None:
        raise ValueError(f"Time token '{time_tok}' from sample '{col}' is not "
                        f"in sample_naming.time_map")
    return time_label, treatment


def sample_metadata(columns, cfg):
    """Build a sample metadata DataFrame from a list of count-matrix columns."""
    rows = []
    for c in columns:
        t, tr = parse_sample(c, cfg)
        rows.append({"sample": c, "time": t, "condition": tr})
    meta = pd.DataFrame(rows).set_index("sample")
    return meta


def time_num(cfg):
    return cfg["sample_naming"]["time_numeric"]


def time_order(cfg):
    return cfg["sample_naming"]["time_order"]


# --------------------------------------------------------------------------
def strip_prefix(series, cfg):
    pref = cfg["inputs"].get("gene_id_prefix_strip", "")
    if pref:
        return series.astype(str).str.replace(f"^{pref}", "", regex=True)
    return series.astype(str)


def drop_configured_samples(df_cols, cfg):
    """Return columns minus any matching drop_samples (exact or prefix)."""
    drop = cfg.get("drop_samples", []) or []
    keep = []
    for c in df_cols:
        if any(c == d or c.startswith(d) for d in drop):
            continue
        keep.append(c)
    return keep


# --------------------------------------------------------------------------
def save_json(obj, path):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=str)


def banner(msg):
    print("\n" + "=" * 72 + f"\n{msg}\n" + "=" * 72, flush=True)
