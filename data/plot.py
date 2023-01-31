import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

data_trotter = np.loadtxt('correlfunc_trotter.csv', delimiter=',',dtype=float)
data_p2tdvp = np.loadtxt('correlfunc_p2tdvp.csv', delimiter=',',dtype=float)
plt.rcParams['text.usetex'] = True
#print(data.shape)
fs=30
fig = plt.figure(figsize=(10,10))
gs  = gridspec.GridSpec(1000, 1000)
#ax2 = plt.subplot(gs[0  :900, 20:920])
ax1 = plt.subplot(gs[0  :600, 20:980])
ax1.plot(data_trotter[:,0],data_trotter[:,1], color="r", label="Trotter", linewidth=4)
ax1.plot(data_p2tdvp[:,0],data_p2tdvp[:,1], "--", color="k", label="p2TDVP", linewidth=4)
ax1.set_ylabel(r"$Re[S^{z}_{L/2}S^{z}_{L/2+1}(t)]$",fontsize=fs)
ax1.set_xlabel(r"$tJ$",fontsize=fs)
ax1.legend(frameon=False, fontsize=fs)
ax1.tick_params(axis="y", labelsize=fs)
ax1.tick_params(axis="x", labelsize=fs)
#ax2.plot(data_trotter[:,0],data_trotter[:,2])
#ax2.plot(data_p2tdvp[:,0],data_p2tdvp[:,2], "--")ax.legend()
plt.savefig("./correl_func.png", format='png', bbox_inches = 'tight')
