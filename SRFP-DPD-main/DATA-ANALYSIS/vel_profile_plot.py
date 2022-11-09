#!/usr/bin/python3
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pylab import cm

# Formatting
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
colors = cm.get_cmap('tab10', 2)

velProfile = np.load("./velProfileData.npy",allow_pickle=True)
print(velProfile.shape)
m = 5
n = 4


    
for j in range(n):
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.plot(np.arange(200)/10.0,np.mean(velProfile[:,j,:], axis = 1),linewidth=2, color = colors(0))
    ax.set_xlabel("Y position",labelpad=10)
    ax.set_ylabel("Avg Velocity (per molecule)",labelpad=10)
    ax.set_xlim(-1, 21)
    ax.set_ylim(-1e-5, 4e-5)
    # ax.set_title("Kinetic energy stability for test "+str(j))
    plt.savefig("velProfile"+str(j+1)+".jpg",dpi=300, transparent=False, bbox_inches='tight')
