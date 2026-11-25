from ase.io import read, write
import os
import random
from ase.build import sort

# 输入 xyz 文件
xyz_file = "mdstruct.xyz"

# 读取所有结构
atoms_list = read(xyz_file, index=":")

n_structures = len(atoms_list)
print(f"Number of structures: {n_structures}")

# 随机选择 20 个结构，并保存其原始编号
n_select = min(20, n_structures)
selected_indices = sorted(random.sample(range(n_structures), n_select))

with open("selected_structures.txt", "w") as f:
    for idx in selected_indices:
        f.write(f"{idx + 1}\n")

print("Selected structures:", [idx + 1 for idx in selected_indices])

# 循环写出并提交随机选中的结构
for idx in selected_indices:
    atoms = atoms_list[idx]
    i = idx + 1
    s = "%03d" % i
    struct = "struct" + s
    print(struct)

    if not os.path.exists(struct):
        os.makedirs(struct)
    else:
        print(f"目录 {struct} 已存在，跳过创建。")
        continue

    s = "cp -f run.in gpumd.sh nep.txt " + struct
    os.system(s)

    os.chdir(struct)
    atoms = sort(atoms)
    write("model.xyz", atoms, format="extxyz")

    s = "bsub < gpumd.sh"
    os.system(s)
    os.chdir("../")
