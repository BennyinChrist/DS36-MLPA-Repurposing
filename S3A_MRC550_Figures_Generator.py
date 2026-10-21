#!/usr/bin/env python3
"""
Supplementary File S3A – MRC-550 Size Standard Figure Generation (Figures 1–5)
"""

import warnings
import zipfile
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

assigned_file = Path("Assigned_Fragments_Output.csv")
if not assigned_file.exists():
    raise FileNotFoundError("Assigned_Fragments_Output.csv not found. Run S1 first.")

df = pd.read_csv(assigned_file).drop_duplicates().reset_index(drop=True)

output_dir = Path("MRC550_Figures")
output_dir.mkdir(exist_ok=True)

rox_all = df[df["Channel"] == "ROX"].copy()
rox_stats_list = []
for nom, g in rox_all.groupby("nominal"):
    rox_stats_list.append({
        "Nominal": nom,
        "Mean_Observed": g["observed"].mean(),
        "Bias": g["observed"].mean() - nom,
        "SD": g["observed"].std(ddof=1) if len(g) > 1 else 0.0
    })
rox_stats = pd.DataFrame(rox_stats_list).sort_values("Nominal").reset_index(drop=True)

retained = [
    58, 68, 78, 87, 98, 107, 119, 123, 126,
    200, 204, 208, 299, 305, 309, 400, 405, 410
]
rox_stats = rox_stats[rox_stats["Nominal"].isin(retained)].copy()

all_mrc = [
    58, 68, 78, 87, 98, 107, 119, 123, 126,
    161, 181, 200, 204, 208, 241, 273, 299, 305, 309,
    374, 400, 405, 410, 440, 458, 460, 480, 500, 506, 511, 550
]

c1_excluded = {161, 181, 241, 273, 374, 440, 458, 550}
c2_excluded = {460, 480}

rox = df[(df["Channel"] == "ROX") & (df["nominal"].isin(retained))].copy()


def plot_electropherogram(data, title, filename):
    fig, ax = plt.subplots(figsize=(12, 5))
    data = data.drop_duplicates(subset=["nominal"]).sort_values("observed").reset_index(drop=True)
    max_height = data["height"].max() if not data.empty else 1000.0

    for i, row in data.iterrows():
        x = float(row["observed"])
        h = float(row["height"])
        nom = float(row["nominal"])

        color = "steelblue" if nom <= 126 else "darkred"
        ax.vlines(x, 0, h, color=color, linewidth=1.5, zorder=2)

        alt_step = (i % 3)
        y_offset = h + (0.05 + alt_step * 0.08) * max_height

        x_nudge = 0.0
        if i > 0:
            prev_x = float(data.iloc[i - 1]["observed"])
            if abs(x - prev_x) < 7.0:
                x_nudge = 3.0 if (i % 2 == 1) else -3.0

        ax.annotate(
            f"{int(nom)} nt",
            xy=(x, h),
            xytext=(x + x_nudge, y_offset),
            textcoords="data",
            ha="center",
            va="bottom",
            fontsize=7,
            fontweight="bold",
            color=color,
            bbox=dict(
                boxstyle="round,pad=0.15",
                facecolor="white",
                edgecolor=color,
                alpha=0.9,
                linewidth=0.6,
            ),
            arrowprops=dict(
                arrowstyle="-", color=color, lw=0.5, shrinkA=0, shrinkB=1
            ),
            zorder=3,
        )

    ax.set_xlabel("Observed Size (nt)", fontsize=10)
    ax.set_ylabel("Peak Height (RFU)", fontsize=10)
    ax.set_title(title, fontsize=11, pad=12)
    ax.set_xlim(40, 460)
    ax.set_ylim(0, max_height * 1.40)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)

    legend_elements = [
        Line2D([0], [0], color="steelblue", lw=2, label="QC zone (58–126 nt, 9/9 fragments)"),
        Line2D([0], [0], color="darkred", lw=2, label="Probe zone (200–410 nt, 9/9 fragments)"),
    ]
    ax.legend(handles=legend_elements, loc="upper right", fontsize=8, frameon=True)

    fig.tight_layout()
    fig.savefig(output_dir / filename, dpi=300, bbox_inches="tight")
    plt.close()


f1_candidates = rox[rox["Well_Type"] == "ROX_Only"]
if not f1_candidates.empty:
    f1_file = f1_candidates["File"].value_counts().idxmax()
    plot_electropherogram(
        rox[rox["File"] == f1_file],
        "Figure 1. Representative electropherogram of the MRC-550 ROX size standard from a ROX-only well",
        "Figure1_ROX_only.png",
    )

f2_candidates = rox[rox["Well_Type"] == "MLPA"]
if not f2_candidates.empty:
    f2_file = f2_candidates["File"].value_counts().idxmax()
    plot_electropherogram(
        rox[rox["File"] == f2_file],
        "Figure 2. Representative electropherogram of the MRC-550 ROX size standard from an MLPA well",
        "Figure2_MLPA.png",
    )

det_rows = []
for nom in all_mrc:
    rox_cnt = df[(df["Channel"] == "ROX") & (df["nominal"] == nom) & (df["Well_Type"] == "ROX_Only")]["File"].nunique()
    mlpa_cnt = df[(df["Channel"] == "ROX") & (df["nominal"] == nom) & (df["Well_Type"] == "MLPA")]["File"].nunique()

    if nom in retained:
        cat = "Retained"
    elif nom in c1_excluded:
        cat = "C1 excluded"
    elif nom in c2_excluded:
        cat = "C2 excluded"
    else:
        cat = "Range cap"

    det_rows.append({"Nominal": nom, "ROX_Only": rox_cnt, "MLPA": mlpa_cnt, "Category": cat})

det = pd.DataFrame(det_rows)
color_map = {
    "Retained": "forestgreen",
    "C1 excluded": "crimson",
    "C2 excluded": "darkorange",
    "Range cap": "purple",
}

fig, ax = plt.subplots(figsize=(13, 5))
x = np.arange(len(det))
width = 0.38

ax.bar(x - width / 2, det["ROX_Only"], width, color=[color_map[c] for c in det["Category"]], edgecolor="black", linewidth=0.3, alpha=0.95)
ax.bar(x + width / 2, det["MLPA"], width, color=[color_map[c] for c in det["Category"]], edgecolor="black", linewidth=0.3, alpha=0.55)

ax.set_xticks(x)
ax.set_xticklabels(det["Nominal"], rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Detection count (Wells)", fontsize=10)
ax.set_xlabel("Nominal size (nt)", fontsize=10)
ax.set_title("Figure 3. Detection counts of MRC-550 fragments across 120 wells", fontsize=11)
ax.axhline(13, color="black", linestyle="--", linewidth=0.9)

legend_handles = [
    Line2D([0], [0], color="forestgreen", lw=4, label="Retained (n=18)"),
    Line2D([0], [0], color="crimson", lw=4, label="C1 excluded (n=8)"),
    Line2D([0], [0], color="darkorange", lw=4, label="C2 excluded (n=2)"),
    Line2D([0], [0], color="purple", lw=4, label="Range cap (n=3)"),
    Line2D([0], [0], color="black", linestyle="--", lw=1, label="C2 threshold (n=13)"),
]
ax.legend(handles=legend_handles, fontsize=8, loc="upper right")
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(output_dir / "Figure3_Detection_Rates.png", dpi=300, bbox_inches="tight")
plt.close()

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
qc_sub = rox_stats[rox_stats["Nominal"] <= 126]
probe_sub = rox_stats[rox_stats["Nominal"] >= 200]

axes[0].scatter(qc_sub["Nominal"], qc_sub["Mean_Observed"], color="steelblue", s=45, label="QC zone (58–126 nt)", zorder=3)
axes[0].scatter(probe_sub["Nominal"], probe_sub["Mean_Observed"], color="darkred", s=45, label="Probe zone (200–410 nt)", zorder=3)

lims = [50, 420]
axes[0].plot(lims, lims, "k--", alpha=0.5, label="Identity (y=x)")
slope, intercept = np.polyfit(rox_stats["Nominal"], rox_stats["Mean_Observed"], 1)
xfit = np.array(lims)
axes[0].plot(xfit, slope * xfit + intercept, color="black", lw=1.2, label=f"Global OLS (slope={slope:.4f})")

r2 = np.corrcoef(rox_stats["Nominal"], rox_stats["Mean_Observed"])[0, 1] ** 2
axes[0].text(0.05, 0.92, f"R² = {r2:.4f}", transform=axes[0].transAxes, fontsize=9, bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

axes[0].set_xlabel("Nominal size (nt)")
axes[0].set_ylabel("Mean observed size (nt)")
axes[0].set_title("Sizing Linearity")
axes[0].set_xlim(lims)
axes[0].set_ylim(lims)
axes[0].legend(fontsize=8, loc="lower right")
axes[0].spines[["top", "right"]].set_visible(False)

colors_bias = ["steelblue" if n <= 126 else "darkred" for n in rox_stats["Nominal"]]
axes[1].bar(rox_stats["Nominal"].astype(str), rox_stats["Bias"], color=colors_bias, edgecolor="black", linewidth=0.4)
mean_global_bias = rox_stats["Bias"].mean()
mean_wr_sd = rox_stats["SD"].mean()

axes[1].axhline(mean_global_bias, color="red", linestyle="-", linewidth=1.2, label=f"Global Mean Bias ({mean_global_bias:+.3f} nt)")
axes[1].axhline(mean_global_bias + mean_wr_sd, color="red", linestyle="--", linewidth=0.9, label=f"±Mean WR_SD ({mean_wr_sd:.3f} nt)")
axes[1].axhline(mean_global_bias - mean_wr_sd, color="red", linestyle="--", linewidth=0.9)

axes[1].set_xlabel("Nominal size (nt)")
axes[1].set_ylabel("Bias (observed − nominal, nt)")
axes[1].set_title("Sizing Bias of Retained Fragments")
axes[1].tick_params(axis="x", rotation=45, labelsize=7)
axes[1].legend(fontsize=8, loc="upper right")
axes[1].spines[["top", "right"]].set_visible(False)

fig.suptitle("Figure 4. Sizing characterization of the 18 retained MRC-550 fragments", fontsize=12, y=1.02)
fig.tight_layout()
fig.savefig(output_dir / "Figure4_Linearity_Bias.png", dpi=300, bbox_inches="tight")
plt.close()

fig, ax = plt.subplots(figsize=(10, 4.5))
ax.bar(rox_stats["Nominal"].astype(str), rox_stats["SD"], color="seagreen", edgecolor="black", linewidth=0.4)
ax.axhline(mean_wr_sd, color="darkred", linestyle=":", linewidth=1.4, label=f"Mean WR_SD = {mean_wr_sd:.3f} nt")

ax.set_xlabel("Nominal size (nt)")
ax.set_ylabel("Within-run SD (nt)")
ax.set_title("Figure 5. Within-fragment sizing precision of the retained MRC-550 peaks", fontsize=11)
ax.tick_params(axis="x", rotation=45, labelsize=7)
ax.legend(fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(output_dir / "Figure5_Precision.png", dpi=300, bbox_inches="tight")
plt.close()

with zipfile.ZipFile("S3A_MRC550_Figures.zip", "w") as zf:
    for png in output_dir.glob("*.png"):
        zf.write(png, png.name)

print(f"S3A Complete: Generated Figures 1–5 in {output_dir}/ and zipped to S3A_MRC550_Figures.zip")
