from itertools import combinations
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import os
from ase.io import write
from pymatgen.io.ase import AseAtomsAdaptor

poscar = "POSCAR"
symprec = 1e-3
save_structures = True

structure = Structure.from_file(poscar)

# O原子编号
o_indices = [
    i for i, site in enumerate(structure)
    if site.species_string == "O"
]

print(f"Total O atoms = {len(o_indices)}")


sga = SpacegroupAnalyzer(structure, symprec=symprec)
symm_ops = sga.get_space_group_operations()

inequivalent_pairs = []

for pair in combinations(o_indices, 2):

    equivalent = False

    frac1 = structure[pair[0]].frac_coords
    frac2 = structure[pair[1]].frac_coords

    for saved_pair in inequivalent_pairs:

        frac3 = structure[saved_pair[0]].frac_coords
        frac4 = structure[saved_pair[1]].frac_coords

        for op in symm_ops:

            p1 = op.operate(frac1) % 1
            p2 = op.operate(frac2) % 1

            cond1 = (
                (abs((p1-frac3+0.5)%1-0.5) < 1e-3).all()
                and
                (abs((p2-frac4+0.5)%1-0.5) < 1e-3).all()
            )

            cond2 = (
                (abs((p1-frac4+0.5)%1-0.5) < 1e-3).all()
                and
                (abs((p2-frac3+0.5)%1-0.5) < 1e-3).all()
            )

            if cond1 or cond2:
                equivalent = True
                break

        if equivalent:
            break

    if not equivalent:
        inequivalent_pairs.append(pair)

# ======================
# 输出结果
# ======================

print()
print("Inequivalent O-vacancy pairs:")
print("-"*50)

for i, pair in enumerate(inequivalent_pairs):

    d = structure.get_distance(pair[0], pair[1])

    print(
        f"{i+1:3d} : "
        f"O({pair[0]})  O({pair[1]})   "
        f"d = {d:.3f} Å"
    )

print()
print(f"Number of inequivalent double O vacancies = {len(inequivalent_pairs)}")


if save_structures:

    all_atoms = []

    for i, pair in enumerate(inequivalent_pairs):

        vac = structure.copy()

        # 倒序删除两个O
        vac.remove_sites(sorted(pair, reverse=True))

        # 保存POSCAR
        vac.to(
            fmt="poscar",
            filename=f"POSCAR_{i+1:03d}"
        )

        # 转为ASE Atoms并保存到列表
        atoms = AseAtomsAdaptor.get_atoms(vac)
        all_atoms.append(atoms)

    # 保存所有结构到一个xyz文件
    write("mdstruct.xyz", all_atoms)

    print(f"Saved {len(inequivalent_pairs)} POSCARs")
    print(f"Saved all structures to mdstruct.xyz")