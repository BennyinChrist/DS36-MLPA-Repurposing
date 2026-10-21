#!/usr/bin/env python3
"""
Supplementary File S2 – Statistical Analysis Module
CLSI EP05-Ed4 variance decomposition and per-fragment performance metrics.
"""

import warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

assigned_file = Path("Assigned_Fragments_Output.csv")
if not assigned_file.exists():
    raise FileNotFoundError("Assigned_Fragments_Output.csv not found. Run S1 first.")

df = pd.read_csv(assigned_file).drop_duplicates().reset_index(drop=True)


def get_run_id(filename):
    fname = Path(filename).name.lower()
    if "095-" in fname:
        return "Run_2"
    elif fname.startswith("mlpa-") or fname.startswith("0") or fname.startswith("1"):
        return "Run_1"
    elif fname.startswith("s"):
        return "Run_3"
    return "Run_1"


df["Run_ID"] = df["File"].apply(get_run_id)

metrics = []
for nom, g in df.groupby("nominal"):
    n_obs = len(g)
    mean_obs = g["observed"].mean()
    bias = mean_obs - nom
    mape = (abs(g["observed"] - nom) / nom * 100).mean()
    wr_sd = g["observed"].std(ddof=1) if n_obs > 1 else 0.0

    rox_g = g[g["Well_Type"] == "ROX_Only"]
    if not rox_g.empty and rox_g["Run_ID"].nunique() > 1:
        run_means = rox_g.groupby("Run_ID")["observed"].mean()
        br_sd = run_means.std(ddof=1)
    else:
        br_sd = 0.0

    cv_pct = (wr_sd / mean_obs * 100) if mean_obs != 0 else 0.0

    metrics.append({
        "Nominal": float(nom),
        "Mean_Observed": round(mean_obs, 3),
        "Bias": round(bias, 3),
        "SD": round(wr_sd, 3),
        "Between_Run_SD": round(br_sd, 3) if not np.isnan(br_sd) else 0.0,
        "MAPE": round(mape, 3),
        "CV_Percent": round(cv_pct, 3),
        "N": int(n_obs),
    })

stats_df = pd.DataFrame(metrics).sort_values("Nominal").reset_index(drop=True)
out_path = Path("Dynamic_Statistical_Metrics.csv")
stats_df.to_csv(out_path, index=False)
print(f"S2 Complete: Saved metrics to {out_path} (Processed {len(stats_df)} nominal fragments)")
