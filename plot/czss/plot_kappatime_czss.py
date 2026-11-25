import numpy as np
import matplotlib.pyplot as plt
from thermo.gpumd.data import load_kappa
from thermo.gpumd.calc import running_ave

# ================= PRB 风格 =================
plt.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 14,
    "axes.linewidth": 1.2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.size": 4,
    "ytick.major.size": 4,
})

def set_ax(ax):
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.minorticks_on()

# ================= 数据 =================
concs = ['0per', '2per', '4per', '6per']
labels = ['0%', '2%', '4%', '6%']

# 平均线颜色
color = '#355C9A'

# 子图编号
letters = ['(a)', '(b)', '(c)', '(d)',
           '(e)', '(f)', '(g)', '(h)']

# ================= 4×2 布局 =================
fig, axes = plt.subplots(4, 2, figsize=(9, 12))

for row, (conc, label) in enumerate(zip(concs, labels)):

    # ================= κxx =================
    ax = axes[row, 0]
    set_ax(ax)

    folders_xx = [f"{conc}/cal{str(i).zfill(3)}" for i in range(1, 6)]

    kxx_list = []

    for folder in folders_xx:

        kappa = load_kappa(folder + "/")

        t = np.arange(1, kappa['kxi'].shape[0] + 1) * 0.001

        kxx = running_ave(kappa['kxi'], t) + running_ave(kappa['kxo'], t)

        kxx_list.append(kxx)

        # 五次独立计算：灰色细线
        ax.plot(t, kxx, color='0.75', linewidth=0.8, alpha=0.9)

    # 平均
    kxx_avg = np.mean(np.array(kxx_list), axis=0)

    ax.plot(t, kxx_avg, color=color, linewidth=2.2, label=label)

    ax.set_xlim(0, 5)
    ax.set_ylim(0, 10)

    if row == 3:
        ax.set_xlabel("Time (ns)")

    ax.set_ylabel(r'$\kappa_{\mathrm{x}}$ (W m$^{-1}$ K$^{-1}$)')

    ax.text(-0.15, 1, letters[2 * row],
            transform=ax.transAxes, fontsize=16)

    ax.text(0.78, 0.85, label,
            transform=ax.transAxes, fontsize=14)

    # ================= κzz =================
    ax2 = axes[row, 1]
    set_ax(ax2)

    folders_zz = [f"{conc}/zz/cal{str(i).zfill(3)}"
                  for i in range(1, 6)]

    kzz_list = []

    for folder in folders_zz:

        kappa = load_kappa(folder + "/")

        kz = running_ave(kappa['kz'], t)

        kzz_list.append(kz)

        # 五次独立计算：灰色细线
        ax2.plot(t, kz, color='0.75', linewidth=0.8, alpha=0.9)

    # 平均
    kzz_avg = np.mean(np.array(kzz_list), axis=0)

    ax2.plot(t, kzz_avg, color=color, linewidth=2.2)

    ax2.set_xlim(0, 5)
    ax2.set_ylim(0, 10)

    if row == 3:
        ax2.set_xlabel("Time (ns)")

    ax2.set_ylabel(r'$\kappa_{\mathrm{z}}$ (W m$^{-1}$ K$^{-1}$)')

    ax2.text(-0.15, 1, letters[2 * row + 1],
             transform=ax2.transAxes, fontsize=16)

    ax2.text(0.78, 0.85, label,
             transform=ax2.transAxes, fontsize=14)

# ================= 布局 =================
plt.tight_layout(h_pad=0.8, w_pad=0.6)

plt.savefig("kappa_8panel.pdf",
            bbox_inches='tight')
plt.savefig("kappa_8panel.png", dpi=400,
            bbox_inches='tight')
plt.show()
