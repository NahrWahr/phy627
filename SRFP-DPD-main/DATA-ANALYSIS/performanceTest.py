import subprocess
import os

m = 5
n = 5

outputEnabled = 0

# numSoParticles Number of particles
# Lx Box length
# Ly Box width
# nt Number of steps
report = open(r"./report.txt", "w")

report.write(" GPU:\n")
report.write("nP init run(50)\n")
subprocess.run(["make", "clean"])
subprocess.run(["make", "NVIDIA=1"])
for j in range(n):
    for i in range(m):
        subprocess.run(["./output/main", str(3000*(j+1)), str(10),
                       str(10*(j+1)), str(50), str(outputEnabled)])
        if outputEnabled:
            os.rename('./outPos.csv', './data/outPosNVIDIA' +
                      str(i)+'_'+str(j)+'.csv')
            os.rename('./outVel.csv', './data/outVelNVIDIA' +
                      str(i)+'_'+str(j)+'.csv')
        lines = open("./time.txt", "r").readlines()
        report.writelines(
            [str(3000*(j+1))," ", lines[0][:-1]," ", lines[1]])
        os.rename('./time.txt', './data/timeNVIDIA'+str(i)+'_'+str(j)+'.csv')

report.write("\n Multicore:\n")
subprocess.run(["make", "clean"])
subprocess.run(["make", "MULTICORE=1"])
for j in range(n):
    for i in range(m):
        subprocess.run(["./output/main", str(3000*(j+1)), str(10),
                       str(10*(j+1)), str(50), str(outputEnabled)])
        if outputEnabled:
            os.rename('./outPos.csv', './data/outPosMulticore' +
                      str(i)+'_'+str(j)+'.csv')
            os.rename('./outVel.csv', './data/outVelMulticore' +
                      str(i)+'_'+str(j)+'.csv')
        lines = open("./time.txt", "r").readlines()
        report.writelines(
            [str(3000*(j+1))," ", lines[0][:-1]," ", lines[1]])
        os.rename('./time.txt', './data/timeMulticore' +
                  str(i)+'_'+str(j)+'.csv')

report.write("\n Serial:\n")
subprocess.run(["make", "clean"])
subprocess.run(["make"])
for j in range(n):
    for i in range(m):
        subprocess.run(["./output/main", str(3000*(j+1)), str(10),
                       str(10*(j+1)), str(50), str(outputEnabled)])
        if outputEnabled:
            os.rename('./outPos.csv', './data/outPosSerial' +
                      str(i)+'_'+str(j)+'.csv')
            os.rename('./outVel.csv', './data/outVelSerial' +
                      str(i)+'_'+str(j)+'.csv')
        lines = open("./time.txt", "r").readlines()
        report.writelines(
            [str(3000*(j+1))," ", lines[0][:-1]," ", lines[1]])
        os.rename('./time.txt', './data/timeSerial'+str(i)+'_'+str(j)+'.csv')
