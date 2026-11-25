from ase.io import read, write
import random

# ===================== 参数 =====================
input_file = "model.xyz"
output_file = "model.xyz"

hg_concentration = 5   # Hg位点中发生Cu/Hg反位的比例，例如10表示10%
# ===============================================

# 读取结构
atoms = read(input_file)

symbols = atoms.get_chemical_symbols()

# 找出Cu和Hg原子的索引
cu_indices = [i for i, s in enumerate(symbols) if s == "Cu"]
hg_indices = [i for i, s in enumerate(symbols) if s == "Hg"]

n_cu = len(cu_indices)
n_hg = len(hg_indices)

print(f"Cu atoms: {n_cu}")
print(f"Hg atoms: {n_hg}")

# 根据Hg浓度计算需要交换的Hg数量
n_swap = round(n_hg * hg_concentration / 100.0)

if n_swap < 1:
    raise ValueError(
        f"当前浓度 {hg_concentration}% 对应的反位数量小于1。"
        f"体系中只有 {n_hg} 个Hg原子。"
    )

if n_swap > n_cu:
    raise ValueError(
        f"需要交换 {n_swap} 个Cu，但体系中只有 {n_cu} 个Cu。"
    )

if n_swap > n_hg:
    raise ValueError(
        f"需要交换 {n_swap} 个Hg，但体系中只有 {n_hg} 个Hg。"
    )

# 随机选择Cu和Hg
selected_cu = random.sample(cu_indices, n_swap)
selected_hg = random.sample(hg_indices, n_swap)

# Cu -> Hg
for i in selected_cu:
    symbols[i] = "Hg"

# Hg -> Cu
for i in selected_hg:
    symbols[i] = "Cu"

atoms.set_chemical_symbols(symbols)

# 保存
write(output_file, atoms)

print()
print(f"Hg antisite concentration: {hg_concentration}%")
print(f"Number of Cu/Hg antisite pairs: {n_swap}")
