import numpy as np
import matplotlib.pyplot as plt

INPUT = "vacancy_distance_summary.dat"
OUTPUT_PNG = "vacancy_distance_distribution.png"
OUTPUT_PDF = "vacancy_distance_distribution.pdf"

# 数据列：
# concentration  mean_dmin_A
# mean_N_lt8  mean_ratio_lt8_percent
# mean_N_8_16  mean_ratio_8_16_percent
# mean_N_gt16  mean_ratio_gt16_percent

data = np.genfromtxt(INPUT, comments="#", dtype=None, encoding="utf-8")

# 兼容只有一行数据的情况
if data.ndim == 0:
    data = np.array([data])

concs = [row[0] for row in data]
ratio_lt8 = np.array([float(row[3]) for row in data])
ratio_8_16 = np.array([float(row[5]) for row in data])
ratio_gt16 = np.array([float(row[7]) for row in data])

plt.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 13,
    "xtick.direction": "in",
    "ytick.direction": "in",
})

fig, ax = plt.subplots(figsize=(6.5, 4.8))

x = np.arange(len(concs))
width = 0.24

bars1 = ax.bar(
    x - width,
    ratio_lt8,
    width,
    label=r"$d_{\mathrm{NN}}<8$ $\AA$"
)

bars2 = ax.bar(
    x,
    ratio_8_16,
    width,
    label=r"$8\leq d_{\mathrm{NN}}\leq16$ $\AA$"
)

bars3 = ax.bar(
    x + width,
    ratio_gt16,
    width,
    label=r"$d_{\mathrm{NN}}>16$ $\AA$"
)

# 在柱子上方标百分比
for bars in [bars1, bars2, bars3]:
    for bar in bars:
        h = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            h,
            f"{h:.1f}%",
            ha="center",
            va="bottom",
        )

ax.set_xlabel("Zn-vacancy concentration (%)")
ax.set_ylabel("Fraction of vacancies (%)")

ax.set_xticks(x)
ax.set_xticklabels([c.replace("%", "") for c in concs])

# 留出百分比文字空间
ymax = max(ratio_lt8.max(), ratio_8_16.max(), ratio_gt16.max())
ax.set_ylim(0, ymax * 1.25)

ax.legend(frameon=False)

plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=600, bbox_inches="tight")
plt.show()
