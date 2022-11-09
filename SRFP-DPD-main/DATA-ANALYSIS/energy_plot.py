#!/usr/bin/python3
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pylab import cm

# Formatting
# mpl.font_manager._rebuild()
mpl.rcParams['font.family'] = 'Coolvetica'
plt.rcParams['font.size'] = 12
plt.rcParams['axes.linewidth'] = 2
colors = cm.get_cmap('tab10', 2)

energy = np.load("./energy.npy",allow_pickle=True)
print(energy.shape)
m = 5
n = 4


    
for j in range(n):
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.semilogx(np.mean(energy[:,j,:], axis = 1)/4800,linewidth=1, color = colors(0))
    ax.set_xlabel("Step (dt = 0.01)",labelpad=10)
    ax.set_ylabel("Kinetic energy (per molecule)",labelpad=10)
    ax.set_ylim(0.4, 1.2)
    # ax.set_title("Kinetic energy stability for test "+str(j))
    plt.savefig("logenergy"+str(j+1)+".jpg",dpi=300, transparent=False, bbox_inches='tight')
