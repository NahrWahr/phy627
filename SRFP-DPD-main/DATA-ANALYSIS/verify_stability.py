from matplotlib import colors
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas.core import frame

Lx = 40
Ly = 20

positions = pd.read_csv('outPos.csv')
ntime, nparticles = positions.shape
data_1 = np.zeros((ntime, nparticles, 2))
for time, row in positions.iterrows():
    for particle, point in enumerate(row):
        x1, x2, x3 = point.split(' ')
        data_1[time, particle, 0] = float(x1[1:])
        data_1[time, particle, 1] = float(x2)

positions = pd.read_csv('outPos_old.csv')
ntime, nparticles = positions.shape
data_2 = np.zeros((ntime, nparticles, 2))
for time, row in positions.iterrows():
    for particle, point in enumerate(row):
        x1, x2, x3 = point.split(' ')
        data_2[time, particle, 0] = float(x1[1:])
        data_2[time, particle, 1] = float(x2)

diff = np.subtract(data_1, data_2)
error = np.sum(np.sqrt(np.multiply(diff, diff)), axis=(1,2))
# print(error.shape)
plt.plot(error)
plt.title("Stability between two simulations with the same parameters.")
plt.ylabel("Mean squared error")
plt.xlabel("Steps")
plt.show()
