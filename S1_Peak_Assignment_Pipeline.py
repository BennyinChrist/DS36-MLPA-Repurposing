#!/usr/bin/env python3
"""
Supplementary File S1 – Peak Assignment Pipeline
Capillary electrophoresis peak assignment and extraction module.
"""

import re
import warnings
import zipfile
from pathlib import Path
import pandas as pd

warnings.filterwarnings("ignore")

# Extract nested zips if present in current working directory
for z in Path(".").glob("*.zip"):
    if any(k in z.name.lower() for k in ["fragment", "labeled", "clean"]):
        with zipfile.ZipFile(z, "r") as zip_ref:
            zip_ref.extractall(".")
        break

files = [f for f in Path(".").rglob("*fragment analysis*.csv") if "raw_data_temp" not in str(f)]
if not files:
    files = [
        f for f in Path(".").rglob("*.csv")
        if "raw_data_temp" not in str(f)
        and not any(x in f.name for x in ["Assigned", "Table", "Script5"])
    ]

if not files:
    raise FileNotFoundError("No raw fragment list CSV files found in directory.")

print(f"Found {len(files)} raw CSV files to process.")

ALL_MRC_550 = [
    58, 68, 78, 87, 98, 107, 119, 123, 126,
    161, 181, 200, 204, 208, 241, 273, 299, 305, 309,
    374, 400, 405, 410, 440, 458, 460, 480, 500, 506, 511, 550
]

RETAINED_MRC_550 = {
    58, 68, 78, 87, 98, 107, 119, 123, 126,
    200, 204, 208, 299, 305, 309, 400, 405, 410
}

P095 = [
    136, 142, 148, 154, 160, 166, 172, 178, 184, 193,
    202, 211, 220, 229, 240, 247, 255, 265, 274, 283,
    291, 301, 310, 319, 337, 346, 355, 364, 382, 391,
    400, 409, 427, 436, 443, 454
]

Y_PROBES = {160, 193, 240, 283}


def infer_well_type(path: Path) -> str:
    pstr = str(path).lower()
    name = path.name.lower()
    if "rox_only" in pstr or "rox-only" in pstr or "rox_only" in name:
        return "ROX_Only"
    if "mlpa" in pstr or "mlpa" in name:
        return "MLPA"
    return "ROX_Only" if name[:3].isdigit() else "MLPA"


def nominal_from_name(name: str):
    if not isinstance(name, str) or name.strip().lower() == "unknown":
        return None
    m = re.search(r"(\d+)\s*nt\b", name, flags=re.IGNORECASE)
    if m:
        return int(m.group(1))
    return None


rows = []
for f in files:
    df_raw = None
    for enc in ["utf-16", "utf-8", "latin1"]:
        try:
            df_raw = pd.read_csv(f, encoding=enc)
            if not df_raw.empty and any(c.lower() in ["length", "height", "size"] for c in df_raw.columns):
                break
        except Exception:
            continue

    if df_raw is None or df_raw.empty:
        continue

    rename_dict = {c: c.capitalize() for c in df_raw.columns if c.lower() in ["length", "height", "name", "size"]}
    df_raw = df_raw.rename(columns=rename_dict)

    if "Size" in df_raw.columns and "Length" not in df_raw.columns:
        df_raw = df_raw.rename(columns={"Size": "Length"})

    if "Name" not in df_raw.columns:
        df_raw["Name"] = "Unknown"

    if "Length" not in df_raw.columns or "Height" not in df_raw.columns:
        continue

    well_type = infer_well_type(f)
    fname = f.name

    work = df_raw[["Name", "Length", "Height"]].dropna(subset=["Length", "Height"]).copy()
    work = work[work["Height"] > 50]
    if work.empty:
        continue

    peaks_sorted = work.sort_values("Length").copy()
    used_indices = set()
    for nom in ALL_MRC_550:
        best_idx, best_dist = None, 6.0
        for idx, r in peaks_sorted.iterrows():
            if idx in used_indices:
                continue
            try:
                dist = abs(float(r["Length"]) - nom)
            except (ValueError, TypeError):
                continue
            if dist < best_dist:
                best_dist = dist
                best_idx = idx
        if best_idx is not None:
            used_indices.add(best_idx)
            r_match = peaks_sorted.loc[best_idx]
            rows.append({
                "observed": float(r_match["Length"]),
                "height": float(r_match["Height"]),
                "nominal": nom,
                "distance": best_dist,
                "Channel": "ROX",
                "File": fname,
                "Well_Type": well_type,
                "Sex_Call": "",
                "Retained": nom in RETAINED_MRC_550,
                "Name": str(r_match["Name"]) if pd.notna(r_match["Name"]) else "Unknown",
                "Assignment": "mrc_nn",
            })

    if well_type == "MLPA":
        for idx, r in work.iterrows():
            name = str(r["Name"]) if pd.notna(r["Name"]) else "Unknown"
            nom = nominal_from_name(name)
            if nom in P095:
                try:
                    obs_val = float(r["Length"])
                    height_val = float(r["Height"])
                except (ValueError, TypeError):
                    continue
                rows.append({
                    "observed": obs_val,
                    "height": height_val,
                    "nominal": nom,
                    "distance": abs(obs_val - nom),
                    "Channel": "FAM",
                    "File": fname,
                    "Well_Type": well_type,
                    "Sex_Call": "",
                    "Retained": False,
                    "Name": name,
                    "Assignment": "named",
                })

columns = ["observed", "height", "nominal", "distance", "Channel", "File", "Well_Type", "Sex_Call", "Retained", "Name", "Assignment"]
result = pd.DataFrame(rows).drop_duplicates().reset_index(drop=True) if rows else pd.DataFrame(columns=columns)

if not result.empty:
    fam_mlpa = result[(result["Channel"] == "FAM") & (result["Well_Type"] == "MLPA")]
    for fname, g in fam_mlpa.groupby("File"):
        y_count = g["nominal"].isin(Y_PROBES).sum()
        sex = "Male" if y_count >= 2 else "Female"
        result.loc[result["File"] == fname, "Sex_Call"] = sex

out_path = Path("Assigned_Fragments_Output.csv")
result.to_csv(out_path, index=False)
print(f"S1 Complete: Successfully processed data. Saved {len(result)} rows to {out_path}")
