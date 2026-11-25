import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.gridspec import GridSpec

# SED 数据路径及前缀
perfect_folder = "perfect"
perfect_prefix = "perfect_czss"
vac_folder = "vac"
vac_prefix = "vac_czss"

# SED 绘图参数
plot_cutoff_freq = 240            # 色散图最高频率 (cm^-1)
plot_interval     = 30            # y 轴刻度间隔 (cm^-1)
plot_color        = "RdBu_r"      # 色散图颜色映射
colorbar_min      = None          # 颜色条最小值 (log10)，None 自动
colorbar_max      = None          # 颜色条最大值 (log10)

lifetime_perfect = "perfect/TOTAL-LORENTZ-Qpoints.Fre_lifetime"
lifetime_vac     = "vac/TOTAL-LORENTZ-Qpoints.Fre_lifetime"

output_filename = "sed.png"
if_show_figures = False

# ========================= 全局绘图样式（PRB 风格） =========================
plt.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 16,
    "axes.linewidth": 1.2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "xtick.minor.size": 3,
    "ytick.minor.size": 3,
    "xtick.major.width": 1.2,
    "ytick.major.width": 1.2,
    "xtick.minor.width": 1.0,
    "ytick.minor.visible": True,
    "legend.frameon": False
})

# ========================= 辅助函数（复用自 plot.py） =========================
def load_sed_data(folder, prefix):
    """从指定文件夹加载 SED 数据、q点距离和高对称点标签"""
    sed_file = os.path.join(folder, prefix + ".SED")
    qpts_file = os.path.join(folder, prefix + ".Qpts")
    freq_file = os.path.join(folder, prefix + ".THz")
    dist_labels_file = os.path.join(folder, prefix + ".Q_distances_and_labels")

    sed_avg = np.loadtxt(sed_file)
    qpoints = np.loadtxt(qpts_file)          # 未使用，保留兼容
    freq_fft = np.loadtxt(freq_file)

    # 读取 q 点距离和高对称点标签
    q_distances = []
    q_labels = {}
    with open(dist_labels_file, "r") as f:
        lines = f.readlines()
        if lines[0].strip() == "Global distances along the paths:":
            q_distances = [float(x) for x in lines[1].split()]
        if len(lines) >= 5 and lines[3].strip() == "High-symmetry points and their distances:":
            for line in lines[4:]:
                parts = line.split(maxsplit=1)
                if len(parts) == 2:
                    distance, label = parts
                    q_labels[float(distance)] = label.strip()
    return sed_avg, q_distances, freq_fft, q_labels

def plot_sed_on_ax(ax, sed_avg, q_distances, freq_fft, q_labels,
                   cutoff_freq_cm1, interval_cm1, cmap):
    """在指定 axes 上绘制单张 SED 色散图（imshow，无 colorbar）"""
    # 对数化
    sed_log = np.log(sed_avg)
    vmin = colorbar_min if colorbar_min is not None else np.floor(sed_log.min())
    vmax = colorbar_max if colorbar_max is not None else np.ceil(sed_log.max())

    # THz -> cm^-1
    freq_cm1 = freq_fft * 33.356
    if cutoff_freq_cm1 is not None:
        cutoff_idx = np.searchsorted(freq_cm1, cutoff_freq_cm1)
        sed_log = sed_log[:cutoff_idx, :]
        freq_cm1 = freq_cm1[:cutoff_idx]

    im = ax.imshow(sed_log, cmap=cmap, interpolation='hanning',
                   aspect='auto', origin='lower', vmin=vmin, vmax=vmax)
    ax.set_xlim(0, len(q_distances) - 1)

    # x 轴：高对称点位置
    pixel_indices = []
    tick_labels = []
    for dist, label in q_labels.items():
        idx = np.argmin(np.abs(np.array(q_distances) - dist))
        pixel_indices.append(idx)
        tick_labels.append(r'$\Gamma$' if label == 'G' else label)
    ax.set_xticks(pixel_indices)
    ax.set_xticklabels(tick_labels)
    ax.get_xticklabels()[0].set_ha('left')

    # y 轴：频率刻度 (cm^-1)
    max_freq = freq_cm1[-1]
    freqs = np.arange(0, np.ceil(max_freq) + 0.01, interval_cm1)
    y_indices = [np.argmin(np.abs(freq_cm1 - f)) for f in freqs]
    ax.set_yticks(y_indices)
    ax.set_yticklabels([f'{f:.0f}' for f in freqs])
    ax.set_ylim(0, y_indices[-1] if cutoff_freq_cm1 is None else cutoff_idx-1)
    ax.set_ylabel('Frequency (cm$^{-1}$)')
    return im

def read_lifetime_data(filepath):
    """读取寿命文件 (频率 THz, 寿命 ps)，返回频率 (cm⁻¹) 和寿命 (ps)"""
    freq_thz = []
    tau_ps = []
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    f_thz = float(parts[0])
                    t_ps = float(parts[1])
                    if t_ps > 0:
                        freq_thz.append(f_thz)
                        tau_ps.append(t_ps)
                except ValueError:
                    continue
    freq_thz = np.array(freq_thz)
    tau_ps = np.array(tau_ps)
    freq_cm = freq_thz * 33.356
    return freq_cm, tau_ps

# ========================= 主程序 =========================
def main():
    # 1. 加载 SED 相关数据
    print("加载完美结构 SED 数据...")
    sed_perf, qdist_perf, freq_perf, labels_perf = load_sed_data(perfect_folder, perfect_prefix)
    print("加载空位结构 SED 数据...")
    sed_vac, qdist_vac, freq_vac, labels_vac = load_sed_data(vac_folder, vac_prefix)

    # 2. 加载寿命数据
    print("加载寿命数据...")
    freq_perf_life, tau_perf = read_lifetime_data(lifetime_perfect)
    freq_vac_life, tau_vac   = read_lifetime_data(lifetime_vac)

    # 固定整图高度为 5 英寸（与原 SED 高度一致）
    fig_height = 5.0
    width_sed1 = 4
    width_sed2 = 2
    width_life = 6    
    total_width = width_sed1 + width_sed2 + width_life

    # 设置 figure 尺寸
    fig = plt.figure(figsize=(total_width, fig_height))
    # 使用 GridSpec 定义三列，宽度比例按实际宽度分配
    gs = GridSpec(1, 3, width_ratios=[width_sed1, width_sed2, width_life],
                  left=0.08, right=0.96, top=0.92, bottom=0.12, wspace=0.3)

    # 创建三个 axes
    ax_sed_perf = fig.add_subplot(gs[0])
    ax_sed_vac  = fig.add_subplot(gs[1], sharey=ax_sed_perf)
    ax_life     = fig.add_subplot(gs[2])

    # 4. 绘制 SED 色散图（完美 + 空位）
    plot_sed_on_ax(ax_sed_perf, sed_perf, qdist_perf, freq_perf, labels_perf,
                   plot_cutoff_freq, plot_interval, plot_color)
    plot_sed_on_ax(ax_sed_vac, sed_vac, qdist_vac, freq_vac, labels_vac,
                   plot_cutoff_freq, plot_interval, plot_color)

    # 右 SED 子图也显示 y 轴刻度（共享 y 轴后默认不显示，这里显式开启）
    ax_sed_vac.tick_params(labelleft=True)

    ax_sed_perf.text(-0.25, 1.04, '(a)', transform=ax_sed_perf.transAxes, fontsize=18, va='top', ha='left')
    ax_sed_vac.text(-0.52, 1.04, '(b)', transform=ax_sed_vac.transAxes, fontsize=18, va='top', ha='left')

    # 5. 绘制寿命散点图
    ax_life.scatter(freq_vac_life, tau_vac, s=24, alpha=0.6, color='#C44E52',
                    label='Vacancy', rasterized=True, linewidths=0)
    ax_life.scatter(freq_perf_life, tau_perf, s=24, alpha=0.6, color='#355C9A',
                    label='Pristine', rasterized=True, linewidths=0)
    ax_life.set_xlabel('Frequency (cm$^{-1}$)')
    ax_life.set_ylabel('Phonon Lifetime (ps)')
    ax_life.set_yscale('log')
    ax_life.set_xlim(0, 250)
    ax_life.set_ylim(0.1, 800)
    ax_life.xaxis.set_major_locator(MultipleLocator(50))
    ax_life.xaxis.set_minor_locator(MultipleLocator(10))
    ax_life.legend(loc='upper right')
    ax_life.text(-0.18, 1.04, '(c)', transform=ax_life.transAxes, fontsize=18, va='top', ha='left')

    # 6. 保存和显示
    plt.savefig(output_filename, dpi=600, bbox_inches='tight')
    if if_show_figures:
        plt.show()
    plt.close(fig)
    print(f"组合图已保存为 {output_filename}")

if __name__ == "__main__":
    main()
