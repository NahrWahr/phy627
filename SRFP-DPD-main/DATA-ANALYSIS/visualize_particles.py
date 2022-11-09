from matplotlib import colors
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from pandas.core import frame

Lx = 40
Ly = 20

data = np.load("/home/debbh/Workspace/SRFP/SRFP-DPD/DATA-ANALYSIS/data/vel40.npy",allow_pickle=True)
# positions = pd.read_csv('./data/outPosNVIDIATest40.csv')
data = data[0]
ntime, nparticles,_ = data.shape
# data = np.zeros((ntime, nparticles, 2))
# for time, row in positions.iterrows():
#     for particle, point in enumerate(row):
#         x1, x2, x3 = point.split(' ')
#         data[time, particle, 0] = float(x1[1:])
#         data[time, particle, 1] = float(x2)
#         # data[time, particle, 2] = float(x3[:-1])

fig = plt.figure(figsize=(40,10))
L = 7.0
ax = plt.axes(xlim=(-1,Lx+1),ylim=(-1,Ly+1))
scatter=ax.scatter(data[0,:,0],data[0,:,1],color="r",s=1)
text = ax.text(0,Ly +1,"",color="r")
ax.set_title("Particle Visalization for "+str(nparticles)+" particles.")
ax.set_xlabel("X pos")
ax.set_ylabel("Y pos")
ax.plot([0,Lx, Lx, 0, 0], [0, 0, Ly, Ly,0],'k')

def step(i):
    temp = np.vstack((data[i,:,0],data[i,:,1]))
    scatter.set_offsets(temp.T)
    text.set_text("t = "+str(i))
    return scatter,

anim = animation.FuncAnimation(fig, step, interval = 20, frames=ntime)

plt.show()
# anim.save("output.mp4")
