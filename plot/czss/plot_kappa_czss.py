import numpy as np
import matplotlib.pyplot as plt
from thermo.gpumd.data import load_kappa
from thermo.gpumd.calc import running_ave

# ================= PRB通用风格 =================
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

# ================= 主图：2x2布局 =================
fig, axes = plt.subplots(2, 2, figsize=(7.5, 6))  # 适当增大尺寸以容纳四子图

# 浓度和颜色（用于时间演化图）
concs = ['0per', '2per', '4per', '6per']
colors = ['#4C72B0', '#55A868', '#DD8452', '#C44E52']

# ========== 第一行：κ 随时间演化 ==========
# (a) κxx vs time
ax = axes[0, 0]
set_ax(ax)
# (b) κzz vs time
ax2 = axes[0, 1]
set_ax(ax2)

t = None  # 时间数组将在第一次读取时确定

for conc, color in zip(concs, colors):
    # ---- κxx ----
    folders_xx = [f"{conc}/cal{str(i).zfill(3)}" for i in range(1, 6)]
    kxx_list = []

    for folder in folders_xx:
        kappa = load_kappa(folder + "/")
        if t is None:
            t = np.arange(1, kappa['kxi'].shape[0] + 1) * 0.001
        kxi_ra = running_ave(kappa['kxi'], t)
        kxo_ra = running_ave(kappa['kxo'], t)
        kxx_list.append(kxi_ra + kxo_ra)

    kxx_avg = np.mean(np.array(kxx_list), axis=0)

    # ---- κzz ----
    folders_zz = [f"{conc}/zz/cal{str(i).zfill(3)}" for i in range(1, 6)]
    kzz_list = []

    for folder in folders_zz:
        kappa = load_kappa(folder + "/")
        kz_ra = running_ave(kappa['kz'], t)
        kzz_list.append(kz_ra)

    kzz_avg = np.mean(np.array(kzz_list), axis=0)

    # 绘图
    axes[0, 0].plot(t, kxx_avg, color=color, linewidth=1.5,
                    label=f"{conc.replace('per', '')}%")
    axes[0, 1].plot(t, kzz_avg, color=color, linewidth=1.5,
                    label=f"{conc.replace('per', '')}%")

# 设置 (a) 子图
ax = axes[0, 0]
ax.set_xlim(0, 5)
ax.set_ylim(0, 8)
ax.set_xlabel('Time (ns)')
ax.set_ylabel(r'$\kappa_{x}$ (W m$^{-1}$ K$^{-1}$)')
ax.text(-0.17, 1, '(a)', transform=ax.transAxes, fontsize=16)
ax.text(0.1, 0.88, 'in-plane', transform=ax.transAxes)

# 设置 (b) 子图
ax = axes[0, 1]
ax.set_xlim(0, 5)
ax.set_ylim(0, 8)
ax.set_xlabel('Time (ns)')
ax.set_ylabel(r'$\kappa_{z}$ (W m$^{-1}$ K$^{-1}$)')
ax.text(-0.17, 1, '(b)', transform=ax.transAxes, fontsize=16)
ax.text(0.1, 0.88, 'out-of-plane', transform=ax.transAxes)

# 图例放在 (b) 中
axes[0, 1].legend(frameon=False, labelspacing=0.2, ncol=2, loc="lower right")

# ========== 第二行：κ 随空位浓度变化 ==========
vac_dirs = ["0per", "2per", "4per", "6per"]
vac_conc = [0, 2, 4, 6]
blue = '#355C9A'

# ---- κxx 浓度依赖 ----
mean_kxx, err_xx = [], []
for vac in vac_dirs:
    folders = [f"{vac}/cal{str(i).zfill(3)}" for i in range(1, 6)]
    vals = []
    for folder in folders:
        kappa = load_kappa(folder + "/")
        t_local = np.arange(1, kappa['kxi'].shape[0] + 1) * 0.001
        kxx = running_ave(kappa['kxi'], t_local) + running_ave(kappa['kxo'], t_local)
        mask = t_local >= (t_local[-1] - 1.0)
        vals.append(np.mean(kxx[mask]))
    arr = np.array(vals)
    mean_kxx.append(arr.mean())
    err_xx.append(np.std(arr, ddof=1))

# ---- κzz 浓度依赖 ----
mean_kzz, err_zz = [], []
for vac in vac_dirs:
    folders = [f"{vac}/zz/cal{str(i).zfill(3)}" for i in range(1, 6)]
    vals = []
    for folder in folders:
        kappa = load_kappa(folder + "/")
        t_local = np.arange(1, kappa['kxi'].shape[0] + 1) * 0.001
        kzz = running_ave(kappa['kz'], t_local)
        mask = t_local >= (t_local[-1] - 1.0)
        vals.append(np.mean(kzz[mask]))
    arr = np.array(vals)
    mean_kzz.append(arr.mean())
    err_zz.append(np.std(arr, ddof=1))

# 子图 (c) κxx vs 浓度
ax = axes[1, 0]
set_ax(ax)
ax.errorbar(vac_conc, mean_kxx, yerr=err_xx, fmt='o', color=blue, markersize=5, capsize=3, linewidth=1.5)
ax.plot(vac_conc, mean_kxx, color=blue, linewidth=1.5)
ax.set_xlabel("Vacancy concentration (%)")
ax.set_ylabel(r"$\kappa_{x}$ (W m$^{-1}$ K$^{-1}$)")
ax.set_xlim(-0.3, 6.3)
ax.set_ylim(0, 8)
ax.text(-0.17, 1, '(c)', transform=ax.transAxes, fontsize=16)

# 子图 (d) κzz vs 浓度
ax = axes[1, 1]
set_ax(ax)
ax.errorbar(vac_conc, mean_kzz, yerr=err_zz, fmt='s', color=blue, markersize=5, capsize=3, linewidth=1.5)
ax.plot(vac_conc, mean_kzz, color=blue, linewidth=1.5)
ax.set_xlabel("Vacancy concentration (%)")
ax.set_ylabel(r"$\kappa_{z}$ (W m$^{-1}$ K$^{-1}$)")
ax.set_xlim(-0.3, 6.3)
ax.set_ylim(0, 8)
ax.text(-0.17, 1, '(d)', transform=ax.transAxes, fontsize=16)

# ================= 保存热导率数据 =================
with open("kappa_data.txt", "w") as f:
    f.write("# Vacancy concentration (%)    kappa_x (W m^-1 K^-1)    Error_x    kappa_z (W m^-1 K^-1)    Error_z\n")

    for c, kx, ex, kz, ez in zip(
        vac_conc, mean_kxx, err_xx, mean_kzz, err_zz
    ):
        f.write(f"{c:>5.1f}    {kx:.6f}    {ex:.6f}    {kz:.6f}    {ez:.6f}\n")

# 整体调整布局
plt.tight_layout(h_pad=0.1, w_pad=0.5)
plt.savefig("kappa.pdf", bbox_inches='tight')
plt.savefig("kappa.png", dpi=600)
plt.show()

