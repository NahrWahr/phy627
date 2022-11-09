#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import gc
import time

startT = 20000
ny = 200
Ly = 20.0
dt = 0.01
binSize = Ly/ny

vel = np.load("./data/vel10.npy", allow_pickle=True)
vel = vel[0]
nt, npar, _ = vel.shape
# print(vel.shape)
m = 4
n = 5
data = np.zeros((nt, m, n))

def getEnergy(i,j):
    vel = np.load("./data/vel"+str(i+1)+str(j)+".npy", allow_pickle=True)
    vel = vel[0]
    energy = np.multiply(vel,vel)
    energy = np.sum(energy, axis=1)
    energy = np.sum(energy, axis=1)
    gc.collect()
    return energy

for i in range (m):
    for j in range(n):
        t_start = time.time()
        energy = getEnergy(i,j)
        t_end = time.time()
        print(i,j,t_end-t_start)
        data[:,i,j] = energy

np.save("energy",data)

