#!/usr/bin/python3
import numpy as np
import pandas as pd
import gc

def fetchData(filepath):
    file = pd.read_csv(filepath)
    ntime, nparticles = file.shape
    data = np.zeros((ntime, nparticles, 2))
    for time, row in file.iterrows():
        for particle, point in enumerate(row):
            x1, x2, x3 = point.split(' ')
            data[time, particle, 0] = float(x1[1:])
            data[time, particle, 1] = float(x2)
    return data, ntime, nparticles
print("Starting data conversion...")
print("0% conversion done",end="",flush=True)
for i in range(1,5):
    for j in range(5):
        data = fetchData("./data/outVelNVIDIATest"+str(i)+str(j)+".csv")
        np.save("vel"+str(i)+str(j),data)
        del data
        gc.collect()
        print("\r"+str((i*5 + j +1)*2)+"% conversion done",end="",flush=True)
        data = fetchData("./data/outPosNVIDIATest"+str(i)+str(j)+".csv")
        np.save("pos"+str(i)+str(j),data)
        del data
        gc.collect()
        print("\r"+str((i*5 + j +1)*2)+"% conversion done",end="",flush=True)

