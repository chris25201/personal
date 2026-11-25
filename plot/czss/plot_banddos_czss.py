import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# ================= PRB风格 =================
plt.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 18,
    "axes.linewidth": 1.2,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True
})

# ================= 读取 DFT band =================
def load_dft_band(filename):
    k, freq = [], []
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            data = list(map(float, line.split()))
            k.append(data[0])
            freq.append(data[1:])
    return np.array(k), np.array(freq)

# ================= 读取 NEP band =================
def load_nep_band(filename):
    data = np.loadtxt(filename)
    k = data[:, 0]
    freq = data[:, 1:]
    return k, freq

# ================= 读取 DOS =================
def load_dos(filename):
    data = np.loadtxt(filename, comments='#')
    freq = data[:, 0]
    total = data[:, 1]
    proj = data[:, 2:]
    return freq, total, proj

# ================= DOS按元素求和 =================
def process_dos_by_element(freq, total, proj):
    Cu = proj[:, 0:2].sum(axis=1)
    Zn = proj[:, 2:3].sum(axis=1)
    Se = proj[:, 3:7].sum(axis=1)
    Sn = proj[:, 7:8].sum(axis=1)
    return freq, total, Cu, Zn, Se, Sn

# ================= 读取数据 =================
k_dft, f_dft = load_dft_band("czss.bands")
k_nep, f_nep = load_nep_band("phonopy_data.txt")

freq_dos, dos_total, dos_proj = load_dos("czss.dos")
freq_dos, dos_total, dos_Cu, dos_Zn, dos_Se, dos_Sn = \
    process_dos_by_element(freq_dos, dos_total, dos_proj)

# NEP 频率转换为 cm⁻¹（如原始单位为 THz）
f_nep *= 33.356

# NEP k-path 归一化并映射到 DFT 的 k 范围
k_nep = (k_nep - k_nep.min()) / (k_nep.max() - k_nep.min())
k_nep = k_nep * (k_dft.max() - k_dft.min()) + k_dft.min()

# ================= 创建上下子图 =================
fig, axes = plt.subplots(2, 1, figsize=(6, 10),
                         gridspec_kw={'height_ratios': [1.2, 1]})

# ================= (a) 声子能带（上） =================
ax_band = axes[0]

# DFT 能带（深蓝）
for i in range(f_dft.shape[1]):
    if i == 0:
        ax_band.plot(k_dft, f_dft[:, i], color='#355C9A', linewidth=2, linestyle='--', label='DFT')
    else:
        ax_band.plot(k_dft, f_dft[:, i], color='#355C9A', linewidth=2, linestyle='--')

# NEP 能带（红）
for i in range(f_nep.shape[1]):
    if i == 0:
        ax_band.plot(k_nep, f_nep[:, i], color='#C44E52', linewidth=2, label='NEP')
    else:
        ax_band.plot(k_nep, f_nep[:, i], color='#C44E52', linewidth=2)

ax_band.set_ylabel("Frequency (cm$^{-1}$)")
ax_band.set_ylim(0, 240)
ax_band.yaxis.set_major_locator(MultipleLocator(30))
ax_band.yaxis.set_minor_locator(MultipleLocator(10))

# 高对称点
kpts = [0.000000, 0.415275, 0.519287, 0.739384, 1.106575, 1.400491,
        1.620588, 1.784923, 2.078567, 2.258609, 2.569873]
labels = ['Γ','X','Y','S','Γ','Z','S','N','P','Y','Z']

ax_band.set_xticks(kpts)
ax_band.set_xticklabels(labels)
ax_band.get_xticklabels()[0].set_ha('left')
ax_band.set_xlim(k_dft.min(), k_dft.max())
ax_band.legend(frameon=False, loc=(0.65,0))
ax_band.text(-0.18, 1, '(a)', transform=ax_band.transAxes, fontsize=20)

# 高对称点竖线
for kp in kpts:
    ax_band.axvline(kp, color='gray', linestyle='--', linewidth=0.5)

# ================= (b) 声子态密度（下） =================
ax_dos = axes[1]

# 总 DOS
ax_dos.plot(freq_dos, dos_total, color='#355C9A', linewidth=2, label='Total')

# 元素投影
ax_dos.plot(freq_dos, dos_Cu, color='#444444', linewidth=1.2, label='Cu')
ax_dos.plot(freq_dos, dos_Zn, color='#d62728', linewidth=1.2, label='Zn')
ax_dos.plot(freq_dos, dos_Sn, color='#17becf', linewidth=1.2, label='Sn')
ax_dos.plot(freq_dos, dos_Se, color='#2ca02c', linewidth=1.2, label='Se')

ax_dos.set_xlabel("Frequency (cm$^{-1}$)")
ax_dos.set_ylabel("PDOS (states/unit cell)")
ax_dos.set_xlim(0, 240)
ax_dos.set_ylim(0, 0.8)   # 自动留出余量
ax_dos.xaxis.set_major_locator(MultipleLocator(30))
ax_dos.xaxis.set_minor_locator(MultipleLocator(10))
ax_dos.yaxis.set_major_locator(MultipleLocator(0.2))
ax_dos.yaxis.set_minor_locator(MultipleLocator(0.04))

ax_dos.legend(frameon=False, loc='upper left')
ax_dos.text(-0.18, 1, '(b)', transform=ax_dos.transAxes, fontsize=20)

# 可选：与上方能带频率范围对齐的灰色参考线（非必须，可删除）
# for y in [0, 60, 120, 180, 240]:
#     ax_dos.axhline(y, color='gray', linestyle=':', linewidth=0.5, alpha=0.5)

# ================= 调整间距 & 输出 =================
#fig.subplots_adjust(hspace=0.1)
plt.tight_layout()
plt.savefig("band_dos.pdf", bbox_inches='tight')
#plt.savefig("band_dos.png", dpi=600)
plt.show()
