from pylab import *
plt.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 12,
    "axes.linewidth": 1.2,
    "xtick.direction": "in",
    "ytick.direction": "in",
})

with open("loss.out", "r") as f:
    lines = f.readlines()
last_line = lines[-1].strip().split()
# 假设格式：step  total_loss  ...   energy_rmse  force_rmse  virial_rmse ...
# 根据说明，第5、6、7列（索引4,5,6）分别是能量、力、维里的 RMSE
energy_rmse = float(last_line[4])   # 0.00087
force_rmse  = float(last_line[5])   # 0.04443
virial_rmse = float(last_line[6])

figure(figsize=(8,7))

subplot(2,2,1)
loss = loadtxt('loss.out')
loglog(loss[:, 1:7])
#loglog(loss[:, 7:9])
xlabel('Generation/100')
ylabel('Loss')
legend(['Total', 'L1', 'L2', 'E-train', 'F-train', 'V-train'], frameon=False, loc='lower left', bbox_to_anchor=(-0.02, -0.02), handlelength=1.0, 
       labelspacing=0.2, handletextpad=0.2, ncol=2, columnspacing=0.8)
text(-0.18, 1.05,'(a)', transform=gca().transAxes, fontsize=14)

subplot(2,2,2)
energy_test = loadtxt('energy_train.out')
plot(energy_test[:, 1], energy_test[:, 0], '.')
plot(linspace(-2.75,-2.3), linspace(-2.75,-2.3), '-')
xlabel('DFT energy (eV/atom)')
ylabel('NEP energy (eV/atom)')
text(-0.18, 1.05,'(b)', transform=gca().transAxes, fontsize=14)
text(0.25, 0.05, f'RMSE = {energy_rmse:.5f} eV/atom',transform=gca().transAxes)

subplot(2,2,3)
force_test = loadtxt('force_train.out')
plot(force_test[:, 3:6], force_test[:, 0:3], '.')
plot(linspace(-12,12), linspace(-12,12), '-')
xlabel('DFT force (eV/A)')
ylabel('NEP force (eV/A)')
legend(['x direction', 'y direction', 'z direction'], frameon=False, handlelength=0.8, handletextpad=0.2)
text(-0.18, 1.05,'(c)', transform=gca().transAxes, fontsize=14)
text(0.25, 0.05, f'RMSE = {force_rmse:.5f} eV/Å',transform=gca().transAxes)

subplot(2,2,4)
virial_test = loadtxt('virial_train.out')
plot(virial_test[:, 6:12], virial_test[:, 0:6], '.')
plot(linspace(-0.8,2.0), linspace(-0.8,2.0), '-')
xlabel('DFT virial (eV/atom)')
ylabel('NEP virial (eV/atom)')
legend(['xx', 'yy', 'zz', 'xy', 'yz', 'zx'], frameon=False, handlelength=0.8, handletextpad=0.2)
text(-0.18, 1.05,'(d)', transform=gca().transAxes, fontsize=14)
text(0.25, 0.05, f'RMSE = {virial_rmse:.5f} eV/atom',transform=gca().transAxes)

tight_layout()
savefig("nep.png", dpi=400)
show()
