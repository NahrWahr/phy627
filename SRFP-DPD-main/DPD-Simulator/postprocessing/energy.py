import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

velocities = pd.read_csv('velocities.csv')
ntime, nparticles = velocities.shape
data = np.zeros((ntime, nparticles, 2))
for time, row in velocities.iterrows():
    for particle, point in enumerate(row):
        x1, x2, x3 = point.split(' ')
        data[time, particle, 0] = float(x1[1:])
        data[time, particle, 1] = float(x2)
        # data[time, particle, 2] = float(x3[:-1])

data = np.multiply(data,data)
energy = np.sum(data,axis=(1,2))
plt.plot(energy)
plt.title("Energy vs Step (KE is not to scale)")
plt.xlabel("step")
plt.ylabel("Kinetic Energy Measure")
plt.show()