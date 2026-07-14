from ase.io import read
import numpy as np
prim = read('model.xyz')
from calorine.calculators import CPUNEP

pot_name = 'nep.txt'
calc = CPUNEP(pot_name)
prim.calc = calc
from calorine.tools import get_force_constants, relax_structure
relax_structure(prim, fmax=0.00001)
phonon = get_force_constants(prim, calc, [9,9,6])
path = 'GXYSGZHNP'
special_points = dict(
    G= [0, 0, 0],
    X= [0, 0, 0.5],
    Y= [-0.130115, 0.130115, 0.5],
    S= [-0.315058, 0.315058, 0.315058],
    Z= [0.5, 0.5, -0.5],
    H= [0.315058, 0.684942, -0.315058],
    N= [0, 0.5, 0],
    P= [0.25, 0.25, 0.25]
)

path_list = []
for start, stop in zip(path[:-1], path[1:]):
    start = special_points[start]
    stop = special_points[stop]
    path_list.append(np.linspace(start, stop, 51))

phonon.run_band_structure(path_list)

band_dict = phonon.get_band_structure_dict()
phonopy_dists = band_dict["distances"]
phonopy_freqs = band_dict["frequencies"]


#phonopy_paths, phonopy_dists, phonopy_freqs, *_ = phonon.get_band_structure()

xticks = [d[0] for d in phonopy_dists] + [phonopy_dists[-1][-1]]


x = np.hstack(phonopy_dists)  # 合并为一维数组
freqs = np.vstack(phonopy_freqs)  # 合并为二维数组 (num_bands, num_points)
# 将每个距离点重复 num_bands 次，并将频率展平
# 将每个距离点重复 num_points 次，并将频率展平
print(freqs.shape)
print(x.shape)
'''x_repeated = np.repeat(x, freqs.shape[1])
freqs_flat = freqs.ravel()  # 直接展平频率数组，不需要转置

# 合并为两列数据
data = np.column_stack((x_repeated, freqs_flat))

# 保存到文件
np.savetxt('phonopy_data.txt', data, fmt='%.6f', header='distance frequency')'''

data = np.column_stack((x, freqs))
with open("phonopy_data.txt", "w") as f:
    f.write("# High-symmetry labels: ")
    f.write(" ".join(list(path)) + "\n")
    f.write("# High-symmetry points: ")
    f.write(" ".join(f"{xx:.6f}" for xx in xticks) + "\n")
    f.write("# 400 rows, 24 columns (each row = one q-point, each col = one phonon branch)\n")
    np.savetxt(f, data, fmt="%.6f")

