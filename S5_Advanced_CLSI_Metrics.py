#!/usr/bin/env python3
"""
Supplementary File S5 – Advanced CLSI-Style Metrics
Total Error EP21 calculation and domain-level linearity assessment.
"""

import warnings
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import linregress

warnings.filterwarnings("ignore")

assigned_file = Path("Assigned_Fragments_Output.csv")
if not assigned_file.exists():
    raise FileNotFoundError("Assigned_Fragments_Output.csv not found. Run S1 first.")

df = pd.read_csv(assigned_file)
df = df.dropna(subset=["observed", "nominal"])
df["nominal"] = df["nominal"].astype(float)
df["observed"] = df["observed"].astype(float)

p095 = [
    136, 142, 148, 154, 160, 166, 172, 178, 184, 193,
    202, 211, 220, 229, 240, 247, 255, 265, 274, 283,
    291, 301, 310, 319, 337, 346, 355, 364, 382, 391,
    400, 409, 427, 436, 443, 454,
]

ep21_rows = []
for nom in p095:
    g = df[(df["Channel"] == "FAM") & (df["nominal"] == nom)]
    if len(g) < 2:
        continue
    biases = g["observed"] - g["nominal"]
    mean_bias = biases.mean()
    sd = biases.std(ddof=1)

    te = abs(mean_bias) + (1.65 * sd)
    mean_obs = g["observed"].mean()
    ep21_rows.append({
        "Nominal_nt": nom,
        "Mean_Observed_nt": round(mean_obs, 3),
        "Bias_nt": round(mean_bias, 3),
        "SD_nt": round(sd, 3),
        "CV_Percent": round(sd / mean_obs * 100, 3) if mean_obs != 0 else np.nan,
        "Total_Error_EP21_nt": round(te, 3),
        "N_Obs": len(g),
    })

ep21 = pd.DataFrame(ep21_rows)
ep21.to_csv("Script5_CLSI_EP21_Probe_Metrics.csv", index=False)


def domain_stats(sub):
    if len(sub) < 3:
        return None
    slope, intercept, r_value, p_value, std_err = linregress(
        sub["nominal"], sub["observed"]
    )
    return {
        "Slope": round(slope, 4),
        "Intercept": round(intercept, 3),
        "R_squared": round(r_value ** 2, 4),
        "P_value": round(p_value, 6),
        "N": len(sub),
    }


lower = df[(df["Channel"] == "ROX") & (df["nominal"] <= 126)]
upper = df[(df["Channel"] == "ROX") & (df["nominal"] >= 200) & (df["nominal"] <= 410)]

lin_rows = []
for name, subset in [("Lower_58_126nt", lower), ("Upper_200_410nt", upper)]:
    res = domain_stats(subset)
    if res:
        res["Domain"] = name
        lin_rows.append(res)

linearity = pd.DataFrame(lin_rows)
linearity.to_csv("Script5_Domain_Linearity_Summary.csv", index=False)
print("S5 Complete: Compiled CLSI EP21 metrics and linearity summary CSVs.")
