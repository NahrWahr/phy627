import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

partions_y = 20.0
Ly = 20.0
n = 1
m = 1

# positions = pd.read_csv('/home/debbh/Workspace/SRFP/SRFP-DPD/DPD-GPU/data/outPosNVIDIATest'+str(1)+str(0)+'.csv')
# ntime, nparticles = positions.shape
# data = np.zeros((m,n,int(partions_y)))
# for j in range(m):
#     for i in range(n):
j = 3
i = 4
positions = pd.read_csv('/home/debbh/Workspace/SRFP/SRFP-DPD/DPD-GPU/data/outPosNVIDIATest'+str(j+1)+str(i)+'.csv')
ntime, nparticles = positions.shape
pos = np.zeros((ntime, nparticles, 2))
for time, row in positions.iterrows():
    for particle, point in enumerate(row):
        x1, x2, x3 = point.split(' ')
        pos[time, particle, 0] = float(x1[1:])
        pos[time, particle, 1] = float(x2)
# print("Extracted Pos data for "+str(j+1)+str(i))

velocities = pd.read_csv('/home/debbh/Workspace/SRFP/SRFP-DPD/DPD-GPU/data/outVelNVIDIATest'+str(j+1)+str(i)+'.csv')
ntime, nparticles = velocities.shape
vel = np.zeros((ntime, nparticles, 2))
for time, row in velocities.iterrows():
    for particle, point in enumerate(row):
        x1, x2, x3 = point.split(' ')
        vel[time, particle, 0] = float(x1[1:])
        vel[time, particle, 1] = float(x2)

# print("Extracted Pos data for "+str(j)+str(i))
velocity_profile = np.zeros(int(partions_y))

for t in range(20000, ntime):
    for p in range(int(partions_y)):
        x1 = np.greater_equal(pos[t, :, 1], p*Ly/partions_y)
        x2 = np.less(pos[t, :, 1], (p+1)*Ly/partions_y)
        particiles_in_partion = np.argwhere(np.logical_and(x1,x2))
        select_velocity = np.squeeze(vel[t, particiles_in_partion, 0])
        velocity_profile[p] += np.sum(select_velocity,axis=0)/ ((ntime-20000) * select_velocity.shape[0])
# print("Calculated Vel profile for for "+str(j)+str(i)+"\n")
plt.plot(velocity_profile)
plt.title("Velocity Profile")
plt.ylabel("Vel")
plt.xlabel("Section")
plt.show()
