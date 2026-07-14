from pylab import *
from ase.io import read, write
from ase.geometry import get_distances
from thermo.gpumd.preproc import add_basis, repeat
import numpy as np


B1 = read('model.xyz')
B1.wrap()

# 找出O原子的索引
o_indices = [i for i, sym in enumerate(B1.get_chemical_symbols()) if sym == 'O']
num_o = len(o_indices)
num_remove = int(0.05 * num_o)  # 空位浓度


positions = B1.get_positions()
removed_indices = []
attempts = 0
max_attempts = 10000

while len(removed_indices) < num_remove and attempts < max_attempts:
    idx = np.random.choice(o_indices)
    if idx in removed_indices:
        attempts += 1
        continue
    removed_indices.append(idx)
    
if len(removed_indices) < num_remove:
    raise RuntimeError(f"无法满足距离条件，仅选择了 {len(removed_indices)} 个 O 空位")

# 从结构中删除这些O原子
mask = np.ones(len(B1), dtype=bool)
mask[removed_indices] = False
B1 = B1[mask]

# 保存最终结构
write("model.xyz", B1)

print(f"总O原子数: {num_o}, 移除: {len(removed_indices)} (5%)")

