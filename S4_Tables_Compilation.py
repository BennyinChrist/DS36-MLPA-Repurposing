#!/usr/bin/env python3
"""
Supplementary File S4 – Dynamic Compilation of Supplementary Tables S2–S5
"""

import warnings
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

assigned_file = Path("Assigned_Fragments_Output.csv")
if not assigned_file.exists():
    raise FileNotFoundError("Assigned_Fragments_Output.csv not found. Run S1 first.")

stats_file = Path("Dynamic_Statistical_Metrics.csv")
if not stats_file.exists():
    raise FileNotFoundError("Dynamic_Statistical_Metrics.csv not found. Run S2 first.")

df = pd.read_csv(assigned_file)
stats = pd.read_csv(stats_file)

retained = [
    58, 68, 78, 87, 98, 107, 119, 123, 126,
    200, 204, 208, 299, 305, 309, 400, 405, 410,
]
s2 = stats[stats["Nominal"].isin(retained)][
    ["Nominal", "Mean_Observed", "Bias", "SD", "N"]
].copy()
s2.columns = ["Nominal_nt", "Mean_Observed_nt", "Bias_nt", "SD_nt", "N"]
s2 = s2.round(3)

p095 = [
    136, 142, 148, 154, 160, 166, 172, 178, 184, 193,
    202, 211, 220, 229, 240, 247, 255, 265, 274, 283,
    291, 301, 310, 319, 337, 346, 355, 364, 382, 391,
    400, 409, 427, 436, 443, 454,
]
s3 = stats[stats["Nominal"].isin(p095)][
    ["Nominal", "Mean_Observed", "Bias", "SD", "MAPE", "N"]
].copy()
s3.columns = ["Nominal_nt", "Mean_Observed_nt", "Bias_nt", "SD_nt", "MAPE_percent", "N"]
s3 = s3.round(3)

well_rows = []
for fname, g in df.groupby("File"):
    rox_g = g[g["Channel"] == "ROX"]
    fam_g = g[g["Channel"] == "FAM"]
    well_type = g["Well_Type"].iloc[0] if not g.empty else ""
    sex = g["Sex_Call"].iloc[0] if "Sex_Call" in g.columns else ""

    n_rox = len(rox_g)
    n_fam = len(fam_g)

    if len(fam_g) >= 2:
        biases = fam_g["observed"] - fam_g["nominal"]
        mean_bias = biases.mean()
        te = abs(mean_bias) + (1.65 * biases.std(ddof=1))
    else:
        mean_bias = np.nan
        te = np.nan

    well_rows.append({
        "File": fname,
        "Well_Type": well_type,
        "Sex_Call": sex,
        "N_ROX": n_rox,
        "N_FAM": n_fam,
        "Mean_Bias_nt": round(mean_bias, 3) if not np.isnan(mean_bias) else "",
        "Total_Error_nt": round(te, 3) if not np.isnan(te) else "",
    })

s4 = pd.DataFrame(well_rows)

overall_bias = (df["observed"] - df["nominal"]).mean()
overall_mape = (np.abs(df["observed"] - df["nominal"]) / df["nominal"] * 100).mean()
mean_sd = stats["SD"].mean()
n_total = len(df)

s5 = pd.DataFrame([
    {"Metric": "Overall Bias (nt)", "Value": round(overall_bias, 3)},
    {"Metric": "Overall MAPE (%)", "Value": round(overall_mape, 3)},
    {"Metric": "Mean within-probe SD (nt)", "Value": round(mean_sd, 3)},
    {"Metric": "Total observations", "Value": n_total},
])

out_xlsx = "Supplementary_Tables_Compiled.xlsx"
try:
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        s2.to_excel(writer, sheet_name="Table_S2", index=False)
        s3.to_excel(writer, sheet_name="Table_S3", index=False)
        s4.to_excel(writer, sheet_name="Table_S4", index=False)
        s5.to_excel(writer, sheet_name="Table_S5", index=False)
except Exception:
    pass

s2.to_csv("Table_S2_Fragment_Shifts.csv", index=False)
s3.to_csv("Table_S3_Probe_Metrics.csv", index=False)
s4.to_csv("Table_S4_Well_Summary.csv", index=False)
s5.to_csv("Table_S5_Validation_Metrics.csv", index=False)

print("S4 Complete: Compiled tables S2–S5 into Excel and CSVs.")
