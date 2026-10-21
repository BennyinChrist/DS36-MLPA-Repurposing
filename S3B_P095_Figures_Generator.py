#!/usr/bin/env python3
"""
Supplementary File S3B – Figure Generation for P095 Probe Characterisation (Figures 6–8)
"""

import warnings
import zipfile
from pathlib import Path
import matplotlib.pyplot as plt
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

p095_noms = [
    136, 142, 148, 154, 160, 166, 172, 178, 184, 193,
    202, 211, 220, 229, 240, 247, 255, 265, 274, 283,
    291, 301, 310, 319, 337, 346, 355, 364, 382, 391,
    400, 409, 427, 436, 443, 454,
]

fam = df[(df["Channel"] == "FAM") & (df["nominal"].isin(p095_noms))].copy()
fam_stats = stats[stats["Nominal"].isin(p095_noms)].copy()

output_dir = Path("P095_Figures")
output_dir.mkdir(exist_ok=True)


def get_domain(n):
    if n < 200:
        return "Short"
    elif n <= 400:
        return "Medium"
    return "Long"


fam_stats["Domain"] = fam_stats["Nominal"].apply(get_domain)

# Figure 6: MAPE Analysis
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
colors = {"Short": "steelblue", "Medium": "seagreen", "Long": "darkred"}

for domain, color in colors.items():
    sub = fam_stats[fam_stats["Domain"] == domain]
    axes[0].scatter(sub["Nominal"], sub["MAPE"], c=color, label=domain, s=45, zorder=3)

overall_mape = fam_stats["MAPE"].mean()
axes[0].axhline(overall_mape, color="gray", linestyle="--", linewidth=1.0, label=f"Overall MAPE = {overall_mape:.3f}%")
axes[0].axhline(1.0, color="black", linestyle=":", linewidth=0.8, label="1% criterion")
axes[0].set_xlabel("Nominal size (nt)")
axes[0].set_ylabel("MAPE (%)")
axes[0].set_title("Per-probe MAPE")
axes[0].legend(fontsize=8)
axes[0].spines[["top", "right"]].set_visible(False)

domain_means = fam_stats.groupby("Domain")["MAPE"].agg(["mean", "std"])
order = ["Short", "Medium", "Long"]
bars = axes[1].bar(
    order,
    domain_means.loc[order, "mean"],
    yerr=domain_means.loc[order, "std"],
    color=[colors[d] for d in order],
    capsize=5,
    edgecolor="black",
    linewidth=0.5,
)
axes[1].axhline(overall_mape, color="gray", linestyle="--", linewidth=1.0)
axes[1].axhline(1.0, color="black", linestyle=":", linewidth=0.8)
axes[1].set_ylabel("MAPE (%)")
axes[1].set_title("Domain-level MAPE")
axes[1].spines[["top", "right"]].set_visible(False)

for bar, dom in zip(bars, order):
    height = bar.get_height()
    axes[1].annotate(
        f"{height:.2f}",
        xy=(bar.get_x() + bar.get_width() / 2, height),
        xytext=(0, 4),
        textcoords="offset points",
        ha="center",
        va="bottom",
        fontsize=8,
    )

fig.suptitle("Figure 6. Mean Absolute Percentage Error (MAPE) of P095 Probes", fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig(output_dir / "Figure6_MAPE.png", dpi=300, bbox_inches="tight")
plt.close()

# Figure 7: Bland-Altman Plot
fam = fam.copy()
fam["diff"] = fam["observed"] - fam["nominal"]
fam["mean_size"] = (fam["observed"] + fam["nominal"]) / 2

mean_bias = fam["diff"].mean()
sd_diff = fam["diff"].std(ddof=1)

loa_upper = mean_bias + (1.96 * sd_diff)
loa_lower = mean_bias - (1.96 * sd_diff)

fig, ax = plt.subplots(figsize=(8, 5.5))
ax.scatter(fam["mean_size"], fam["diff"], alpha=0.45, s=14, c="steelblue", zorder=2)
ax.axhline(mean_bias, color="darkred", linewidth=1.4, label=f"Mean bias = {mean_bias:+.3f} nt")
ax.axhline(loa_upper, color="gray", linestyle="--", linewidth=1.0, label=f"LoA upper = {loa_upper:+.3f} nt")
ax.axhline(loa_lower, color="gray", linestyle="--", linewidth=1.0, label=f"LoA lower = {loa_lower:+.3f} nt")
ax.set_xlabel("Mean size (nt)")
ax.set_ylabel("Observed − Nominal (nt)")
ax.set_title("Figure 7. Bland–Altman Plot – P095 Probe Sizing")
ax.legend(fontsize=8, loc="upper right")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(output_dir / "Figure7_BlandAltman.png", dpi=300, bbox_inches="tight")
plt.close()

# Figure 8: Background Noise Floor Histogram
rox_only = df[(df["Well_Type"] == "ROX_Only") & (df["height"] < 500)]
heights = rox_only["height"].values if not rox_only.empty else np.array([])

fig, ax = plt.subplots(figsize=(8, 4.5))
if len(heights) > 0:
    ax.hist(heights, bins=40, color="steelblue", edgecolor="white", alpha=0.85)
    mean_h = float(np.mean(heights))
    ax.axvline(
        mean_h,
        color="darkred",
        linestyle="--",
        linewidth=1.3,
        label=f"Mean = {mean_h:.1f} RFU",
    )
ax.set_xlabel("Peak height (RFU)")
ax.set_ylabel("Count")
ax.set_title("Figure 8. Background peaks (ROX-only wells, height < 500 RFU)")
ax.legend(fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(output_dir / "Figure8_Background.png", dpi=300, bbox_inches="tight")
plt.close()

with zipfile.ZipFile("S3B_P095_Figures.zip", "w") as zf:
    for png in output_dir.glob("*.png"):
        zf.write(png, png.name)

print(f"S3B Complete: Generated Figures 6–8 in {output_dir}/ and zipped to S3B_P095_Figures.zip")
