import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation


Lx = 40
Ly = 20

positions = pd.read_csv('./data/outPosNVIDIATest40.csv')
ntime, nparticles = positions.shape
xpos = np.zeros((ntime, nparticles))
ypos = np.zeros((ntime, nparticles))

xedges = np.arange(0,Lx+1, 1)
yedges = np.arange(0,Ly+1, 1)


for time, row in positions.iterrows():
    for particle, point in enumerate(row):
        x1, x2, x3 = point.split(' ')
        xpos[time, particle] = float(x1[1:])
        ypos[time, particle] = float(x2)

print("Data Loaded")

density = np.zeros((ntime, Lx, Ly))
for time in range(ntime):
    density[time,:,:],_,_ = np.histogram2d(xpos[time,:], ypos[time,:], bins=(xedges,yedges)) 

print("Density ")

fig = plt.figure(figsize=(40,20))

ax = plt.axes()
img=ax.imshow(density[0,:,:].T)
# text = ax.text(0,Ly +1,"",color="r")
# ax.set_title("Particle Visalization for "+str(nparticles)+" particles.")
# ax.set_xlabel("X pos")
# ax.set_ylabel("Y pos")
# ax.plot([0,Lx, Lx, 0, 0], [0, 0, Ly, Ly,0],'k')

def step(i):
    # temp = np.vstack((data[i,:,0],data[i,:,1]))
    img.set_array(np.squeeze(density[i,:,:].T))
    # text.set_text("t = "+str(i))
    return img,

anim = animation.FuncAnimation(fig, step, interval = 20, frames=ntime)

plt.show()