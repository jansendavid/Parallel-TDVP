import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('overlap_trotter.csv', delimiter=',',dtype=float)
print(data.shape)
fig = plt.figure()
gs  = gridspec.GridSpec(1000, 1000)
ax2 = plt.subplot(gs[0  :445, 20:920])
ax1 = plt.subplot(gs[544  :990, 20:920])
ax1.plot(data)
