import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

partions_y = 20
Ly = 10

positions = pd.read_csv('positions.csv')
ntime, nparticles = positions.shape
pos = np.zeros((ntime, nparticles, 2))
for time, row in positions.iterrows():
    for particle, point in enumerate(row):
        x1, x2, x3 = point.split(' ')
        pos[time, particle, 0] = float(x1[1:])
        pos[time, particle, 1] = float(x2)

velocities = pd.read_csv('velocities.csv')
ntime, nparticles = velocities.shape
vel = np.zeros((ntime, nparticles, 2))
for time, row in velocities.iterrows():
    for particle, point in enumerate(row):
        x1, x2, x3 = point.split(' ')
        vel[time, particle, 0] = float(x1[1:])
        vel[time, particle, 1] = float(x2)

velocity_profile = np.zeros((partions_y, 2))

for t in range(500, ntime):
    for p in range(partions_y):
        particiles_in_partion = np.argwhere(np.logical_and(pos[t, :, 0] >= p*Ly/partions_y, pos[t, :, 0] < (p+1)*Ly/partions_y))
        select_velocity = np.squeeze(vel[t, particiles_in_partion, :])
        velocity_profile[p, :] += np.sum(select_velocity,axis=0)/ ((ntime-500) * select_velocity.shape[0])

plt.plot(velocity_profile[:,0])
plt.title("Velocity Profile")
plt.ylabel("Vel")
plt.xlabel("Section")
plt.show()
