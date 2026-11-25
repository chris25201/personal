from ase.io import read, write
import os
from ase.build import sort

# 输入 xyz 文件
xyz_file = "make_dataset.xyz"

# 读取所有结构
atoms_list = read(xyz_file, index=":")

# 循环写出 POSCAR
for i, atoms in enumerate(atoms_list, start=1):
    atoms = sort(atoms)
    write(f"POSCAR_{i:03d}", atoms, format="vasp")
#    print(f"已生成 {filename}")

