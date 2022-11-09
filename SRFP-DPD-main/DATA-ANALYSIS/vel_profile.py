#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import gc

startT = 20000
ny = 200
Ly = 20.0
dt = 0.01
binSize = Ly/ny

m = 4
n = 5

def getVelProfile(i,j):
    vel = np.load("./data/vel"+str(i+1)+str(j)+".npy", allow_pickle=True)
    pos = np.load("./data/pos"+str(i+1)+str(j)+".npy", allow_pickle=True)

    vel = vel[0]
    pos = pos[0]

    gc.collect()

    nt, npar, _ = vel.shape
    

    binnedVelocity = np.zeros((nt, ny))

    for t in range(startT, nt):
        for y in range(ny):
            validPos = np.logical_and(
                y*binSize <= pos[t, :, 1], pos[t, :, 1] < (y+1)*binSize)
            binnedVelocity[t, y] = np.sum(vel[t, validPos, 0])

    velProfile = np.sum(binnedVelocity, axis=0) * dt / ((nt - startT)*npar)
    return velProfile


data = np.zeros((ny, m, n))
print("0% done",end="", flush=True)
for i in range(m):
    for j in range(n):
        data[:,i,j] = getVelProfile(i,j)
        print("\r"+str((i*m + j+1)*100/20.0)+"% done",end="",flush=True)
np.save("velProfileData",data)
