#!/usr/bin/env python3
"""
Visualize IDR, D-rich segment, and position-specific analysis results.
"""

import matplotlib.pyplot as plt
import numpy as np

# results from analyze_idr_drich_position.py
# (hardcoded from output to avoid re-running)

BIN_NAMES = ["N-term", "Early-mid", "Middle", "Late-mid", "C-term"]

# D+E (acidic) by position
rbd_de = [10.29, 11.92, 11.85, 13.10, 12.28]
nonrbd_de = [12.08, 11.17, 11.05, 11.16, 10.02]

# K+R (basic) by position
rbd_kr = [11.41, 11.79, 12.92, 10.27, 13.81]
nonrbd_kr = [14.04, 12.29, 12.58, 10.03, 14.40]

# net charge by position
rbd_net = [1.12, -0.13, 1.07, -2.83, 1.53]
nonrbd_net = [1.96, 1.12, 1.53, -1.13, 4.38]

# summary stats
summary_data = {
    "D/E segment presence": {"RBD": 51.2, "non-RBD": 35.7},
    "Disorder bias (score > 0)": {"RBD": 71.4, "non-RBD": 79.5},
    "Mean disorder score": {"RBD": 0.0663 * 100, "non-RBD": 0.1070 * 100},  # scaled for visibility
}

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle("RBD vs non-RBD: Position-Specific Composition Analysis", fontsize=14, fontweight="bold")

# plot 1: D+E by position
ax1 = axes[0, 0]
x = np.arange(len(BIN_NAMES))
width = 0.35
ax1.bar(x - width/2, rbd_de, width, label="RBD", color="#E74C3C", alpha=0.8)
ax1.bar(x + width/2, nonrbd_de, width, label="non-RBD", color="#3498DB", alpha=0.8)
ax1.set_ylabel("Frequency (%)")
ax1.set_title("Acidic (D+E) Residues by Position")
ax1.set_xticks(x)
ax1.set_xticklabels(BIN_NAMES, rotation=15, ha="right")
ax1.legend()
ax1.axhline(y=np.mean(rbd_de), color="#E74C3C", linestyle="--", alpha=0.5)
ax1.axhline(y=np.mean(nonrbd_de), color="#3498DB", linestyle="--", alpha=0.5)

# annotate key differences
for i, (r, n) in enumerate(zip(rbd_de, nonrbd_de)):
    if abs(r - n) > 1.0:
        ax1.annotate(f"{r-n:+.1f}%", (i, max(r, n) + 0.3), ha="center", fontsize=8, color="red")

# plot 2: K+R by position
ax2 = axes[0, 1]
ax2.bar(x - width/2, rbd_kr, width, label="RBD", color="#E74C3C", alpha=0.8)
ax2.bar(x + width/2, nonrbd_kr, width, label="non-RBD", color="#3498DB", alpha=0.8)
ax2.set_ylabel("Frequency (%)")
ax2.set_title("Basic (K+R) Residues by Position")
ax2.set_xticks(x)
ax2.set_xticklabels(BIN_NAMES, rotation=15, ha="right")
ax2.legend()

# annotate key differences
for i, (r, n) in enumerate(zip(rbd_kr, nonrbd_kr)):
    if abs(r - n) > 1.0:
        ax2.annotate(f"{r-n:+.1f}%", (i, max(r, n) + 0.3), ha="center", fontsize=8, color="red")

# plot 3: net charge by position
ax3 = axes[1, 0]
ax3.plot(BIN_NAMES, rbd_net, "o-", color="#E74C3C", label="RBD", linewidth=2, markersize=8)
ax3.plot(BIN_NAMES, nonrbd_net, "s-", color="#3498DB", label="non-RBD", linewidth=2, markersize=8)
ax3.axhline(y=0, color="gray", linestyle="--", alpha=0.5)
ax3.set_ylabel("Net Charge (K+R - D-E) %")
ax3.set_title("Net Charge Profile Along Domain")
ax3.set_xticks(range(len(BIN_NAMES)))
ax3.set_xticklabels(BIN_NAMES, rotation=15, ha="right")
ax3.legend()
ax3.fill_between(range(len(BIN_NAMES)), rbd_net, nonrbd_net, alpha=0.2, color="gray")

# highlight the late-mid negative region
ax3.annotate("More acidic\nin RBDs", (3, -2.0), fontsize=9, ha="center",
             arrowprops=dict(arrowstyle="->", color="red"), xytext=(3.5, 0.5))

# plot 4: summary comparison
ax4 = axes[1, 1]
categories = list(summary_data.keys())
rbd_vals = [summary_data[k]["RBD"] for k in categories]
nonrbd_vals = [summary_data[k]["non-RBD"] for k in categories]

x4 = np.arange(len(categories))
ax4.bar(x4 - width/2, rbd_vals, width, label="RBD", color="#E74C3C", alpha=0.8)
ax4.bar(x4 + width/2, nonrbd_vals, width, label="non-RBD", color="#3498DB", alpha=0.8)
ax4.set_ylabel("Percentage")
ax4.set_title("Summary: D/E Segments & Disorder")
ax4.set_xticks(x4)
ax4.set_xticklabels(categories, rotation=15, ha="right", fontsize=9)
ax4.legend()

# add text annotation for key finding
ax4.annotate("51% of RBDs have\nD/E-rich segments\nvs 36% non-RBD",
             (0, 55), fontsize=9, ha="center", style="italic",
             bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.5))

plt.tight_layout()
plt.savefig("/storage/kiran-stuff/protein-semantic-collapse/results/position_analysis.png", dpi=150)
plt.savefig("/storage/kiran-stuff/protein-semantic-collapse/results/position_analysis.pdf")
print("Saved: results/position_analysis.png and .pdf")

# ========================================
# FIGURE 2: literature comparison
# ========================================
fig2, ax = plt.subplots(figsize=(10, 6))

# our findings vs literature expectations
findings = [
    ("D enrichment in RBDs", "+0.98%", "✓", "Supports D/E-rich IDR hypothesis"),
    ("K/R at interface", "Expected high", "~", "Domain-wide shows slight depletion"),
    ("D/E segments", "51% vs 36%", "✓", "RBDs have more autoinhibitory regions"),
    ("Late-mid acidic", "+1.94%", "✓", "Autoinhibitory elements near C-term"),
    ("Disorder score", "Lower in RBD", "✗", "Opposite: RBDs more ordered"),
]

# create table
ax.axis("off")
table = ax.table(
    cellText=[(f[0], f[1], f[2], f[3]) for f in findings],
    colLabels=["Finding", "Value", "Match", "Interpretation"],
    loc="center",
    cellLoc="left",
    colWidths=[0.25, 0.15, 0.08, 0.52]
)
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.8)

# color cells based on match
for i, f in enumerate(findings):
    if f[2] == "✓":
        table[(i+1, 2)].set_facecolor("#90EE90")  # light green
    elif f[2] == "✗":
        table[(i+1, 2)].set_facecolor("#FFB6C1")  # light red
    else:
        table[(i+1, 2)].set_facecolor("#FFFACD")  # light yellow

ax.set_title("Our Data vs Literature Predictions (Wang et al. 2025)", fontsize=12, fontweight="bold", pad=20)

plt.tight_layout()
plt.savefig("/storage/kiran-stuff/protein-semantic-collapse/results/literature_comparison.png", dpi=150)
print("Saved: results/literature_comparison.png")

print("\nKey findings summary:")
print("=" * 60)
print("1. D/E-rich segments: RBDs have MORE (51% vs 36%)")
print("   → SUPPORTS autoinhibition hypothesis")
print()
print("2. Position bias: D+E enriched at late-mid (+1.94%)")
print("   and C-term (+2.26%) in RBDs")
print("   → Suggests autoinhibitory elements near C-terminus")
print()
print("3. Net charge: RBDs more negative, especially late-mid")
print("   (−2.83 vs −1.13) and C-term (+1.53 vs +4.38)")
print("   → Acidic patches that could compete with RNA")
print()
print("4. UNEXPECTED: RBDs have LOWER disorder scores")
print("   → Domain cores are more ordered, IDRs may be")
print("   in flanking regions outside our domain boundaries")
