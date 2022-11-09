#!/usr/bin/python3
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pylab import cm

# Formatting
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
colors = cm.get_cmap('tab10', 2)

density = np.load("./densityProfile.npy",allow_pickle=True)
print(density.shape)
m = 5
n = 4
for j in range(n):
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # print(density.shape)
    densityProfile = np.mean(np.mean(density[20000:,j,:,:],axis=0),axis=0) +0.3
    densityProfile[-2:] = None
    ax.plot(np.arange(200)/10.0,densityProfile)
    ax.set_xlabel("Y position",labelpad=10)
    ax.set_ylabel("Density",labelpad=10)
    plt.savefig("densityProfile"+str(j+1)+".jpg",dpi=300, transparent=False, bbox_inches='tight')

