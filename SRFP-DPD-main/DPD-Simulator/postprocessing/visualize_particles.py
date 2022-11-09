from matplotlib import colors
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from pandas.core import frame

positions = pd.read_csv('positions.csv')
ntime, nparticles = positions.shape
data = np.zeros((ntime, nparticles, 2))
for time, row in positions.iterrows():
    for particle, point in enumerate(row):
        x1, x2, x3 = point.split(' ')
        data[time, particle, 0] = float(x1[1:])
        data[time, particle, 1] = float(x2)
        # data[time, particle, 2] = float(x3[:-1])

fig = plt.figure(figsize=(40,10))
L = 7.0
ax = plt.axes(xlim=(0,40),ylim=(0,10))
scatter=ax.scatter(data[0,:,0],data[0,:,1])
text = ax.text(0,10.1,"",color="r")
ax.set_title("Particle Visalization for "+str(nparticles)+" particles.")
ax.set_xlabel("X pos")
ax.set_ylabel("Y pos")

def step(i):
    temp = np.vstack((data[i,:,0],data[i,:,1]))
    scatter.set_offsets(temp.T)
    text.set_text("t = "+str(i))
    return scatter,

anim = animation.FuncAnimation(fig, step, interval = 20, frames=ntime)

plt.show()
# anim.save("output.mp4")
