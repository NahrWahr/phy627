import collections
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

gpudata_init       = np.zeros(5)
multicoredata_init = np.zeros(5)
serialdata_init    = np.zeros(5)
gpudata_run        = np.zeros(5)
multicoredata_run  = np.zeros(5)
serialdata_run     = np.zeros(5)


df   = pd.read_csv('./perfReport.txt', sep=' ')
data = np.array(df)
for i in range(5): 
    gpudata_init       [i] = np.mean(data[5*i:5*(i+1), 1])
    multicoredata_init[i]  = np.mean(data[5*i+25:5*(i+1)+25, 1])
    serialdata_init    [i] = np.mean(data[5*i+50:5*(i+1)+50, 1])
    gpudata_run        [i] = np.mean(data[5*i:5*(i+1), 2])
    multicoredata_run  [i] = np.mean(data[5*i+25:5*(i+1)+25, 2])
    serialdata_run     [i] = np.mean(data[5*i+50:5*(i+1)+50, 2])

x     = np.arange(5)
width = 0.2

# plt.bar(x-width, serialdata_init, width, color="grey")
# plt.bar(x,    multicoredata_init, width, color="orange")
# plt.bar(x+width,    gpudata_init, width, color="green")
# plt.xticks(x, [3000*(i+1) for i in range(5)])
# plt.xlabel("Number of particles")
# plt.ylabel("Time taken in ms (Avg over 5 trials)")
# # plt.title("Time taken for initalization")
# plt.legend(["Single-thread", "Multi-threaded", "GPU Acclerated"])
# plt.show()
plt.figure(figsize=(8,6))
plt.bar(x-width, serialdata_run, width, color="grey")
plt.bar(x,    multicoredata_run, width, color="orange")
plt.bar(x+width,    gpudata_run, width, color="green")
plt.xticks(x, [3000*(i+1) for i in range(5)])
plt.xlabel("Number of particles (50 step iteration)")
plt.ylabel("Time taken in ms (Avg over 5 trials)")
plt.title("Time taken to run")
plt.legend(["Single-thread", "Multi-threaded", "GPU Acclerated"])
plt.savefig("perfReport.jpg",dpi=300)
plt.show()
