#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import gc
import time

startT = 20000
bins = 10
Ly = 20.0
Lx = 40.0
dt = 0.01
binSize = Ly/bins

vel = np.load("./data/pos10.npy", allow_pickle=True)
vel = vel[0]
nt, npar, _ = vel.shape
# print(vel.shape)
m = 4
n = 5
data = np.zeros((nt, m, n, int(bins*Ly)))
def getDensityProfile(i,j):
    pos = np.load("./data/pos"+str(i+1)+str(j)+".npy", allow_pickle=True)
    pos = pos[0]
    for t in range(nt):
        for p in range(npar):
            y = int(np.floor(pos[t,p,1]*bins))
            # print(y)
            data[t,i,j,y-1] += 1
print("0% done",end="",flush=True)
for i in range(m):
    for j in range(n):
        start = time.time()
        getDensityProfile(i,j)
        end = time.time()
        print("\r"+str((i*5+j+1)*100.0/20)+"% done in time "+str(end-start),end="",flush=True)

np.save("densityProfile",data)

