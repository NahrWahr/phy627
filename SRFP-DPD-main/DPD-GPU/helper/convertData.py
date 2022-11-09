import numpy as np;
import pandas as pd

print("Loading Pos Data...")

positions = pd.read_csv('/home/debbh/Workspace/SRFP/SRFP-DPD/DPD-GPU/data/outPosNVIDIATest40.csv')
ntime, nparticles = positions.shape
pos = np.zeros((ntime, nparticles, 2))
for time, row in positions.iterrows():
    for particle, point in enumerate(row):
        x1, x2, x3 = point.split(' ')
        pos[time, particle, 0] = float(x1[1:])
        pos[time, particle, 1] = float(x2)

print("\rLoading Vel Data...",flush=True)

velocities = pd.read_csv('/home/debbh/Workspace/SRFP/SRFP-DPD/DPD-GPU/data/outVelNVIDIATest40.csv')
ntime, nparticles = velocities.shape
vel = np.zeros((ntime, nparticles, 2))
for time, row in velocities.iterrows():
    for particle, point in enumerate(row):
        x1, x2, x3 = point.split(' ')
        vel[time, particle, 0] = float(x1[1:])
        vel[time, particle, 1] = float(x2)

print("\rSaving Data into files...",flush=True)

np.savez_compressed('./data/T10', pos=pos, vel=vel)

print("\rData converstion and save  successful!         ",flush=True)
