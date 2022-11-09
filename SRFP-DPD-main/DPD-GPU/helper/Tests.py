import subprocess
import os

subprocess.run(["make", "clean"])
subprocess.run(["make", "NVIDIA=1"])
# subprocess.run(["make", "SERIAL=1"])

n = 5

'''
# Test 1
numSoParticles = str(3200)
Lx             = str(40)
Ly             = str(20)
Lz             = str(0)
Fe             = str(0.0)
nt             = str(10000)
outputEnabled  = str(1)
isSlipBoundary = str(0)
isPeriodic     = str(0)

### Demonstrate: 
* Energy stabilizes over time.
* Particle density stabilizes over time.
'''


numSoParticles = str(2400)
Lx             = str(40)
Ly             = str(20)
Lz             = str(0)
Fe             = str(0.0)
xVel           = str(0.0)
nt             = str(40000)
outputEnabled  = str(1)
isSlipBoundary = str(1)
isPeriodic     = str(0)

for i in range(n): 
    subprocess.run(["./output/main", numSoParticles, Lx, Ly, Lz, Fe, xVel,
                   nt, outputEnabled, isSlipBoundary, isPeriodic, str(2020 + i)])
    if bool(outputEnabled): 
        os.rename('./outPos.csv', './data/outPosNVIDIATest1' +
                  str(i)+'.csv')
        os.rename('./outVel.csv', './data/outVelNVIDIATest1' +
                  str(i)+'.csv')

'''
# Test 2
numSoParticles = str(3200)
Lx             = str(40)
Ly             = str(20)
Lz             = str(0)
Fe             = str(0.0)
xVel           = str(0.5)
nt             = str(10000)
outputEnabled  = str(1)
isSlipBoundary = str(1)
isPeriodic     = str(1)
Velocity is boosted in x direction

### Demonstrate:
* Energy stabilizes over time.
* Particle density stabilizes over time.
'''

numSoParticles = str(2400)
Lx             = str(40)
Ly             = str(20)
Lz             = str(0)
Fe             = str(0.0)
xVel           = str(0.5)
nt             = str(40000)
outputEnabled  = str(1)
isSlipBoundary = str(1)
isPeriodic     = str(1)

for i in range(n): 
    subprocess.run(["./output/main", numSoParticles, Lx, Ly, Lz, Fe, xVel,
                   nt, outputEnabled, isSlipBoundary, isPeriodic, str(2020 + i)])
    if bool(outputEnabled): 
        os.rename('./outPos.csv', './data/outPosNVIDIATest2' +
                  str(i)+'.csv')
        os.rename('./outVel.csv', './data/outVelNVIDIATest2' +
                  str(i)+'.csv')

'''
# Test 3
numSoParticles = str(3200)
Lx             = str(40)
Ly             = str(20)
Lz             = str(0)
Fe             = str(0.0)
xVel           = str(0.5)
nt             = str(10000)
outputEnabled  = str(1)
isSlipBoundary = str(0)
isPeriodic     = str(1)
Velocity is boosted in x direction

### Demonstrate:
* Energy decreases over time.
* Particle density stabilizes over long time.
* Velocity profile after long time.
'''


numSoParticles = str(2400)
Lx             = str(40)
Ly             = str(20)
Lz             = str(0)
Fe             = str(0.0)
xVel           = str(0.5)
nt             = str(40000)
outputEnabled  = str(1)
isSlipBoundary = str(0)
isPeriodic     = str(1)

for i in range(n): 
    subprocess.run(["./output/main", numSoParticles, Lx, Ly, Lz, Fe, xVel,
                   nt, outputEnabled, isSlipBoundary, isPeriodic, str(2020 + i)])
    if bool(outputEnabled): 
        os.rename('./outPos.csv', './data/outPosNVIDIATest3' +
                  str(i)+'.csv')
        os.rename('./outVel.csv', './data/outVelNVIDIATest3' +
                  str(i)+'.csv')

'''
# Test 4
numSoParticles = str(3200)
Lx             = str(40)
Ly             = str(20)
Lz             = str(0)
Fe             = str(0.28)
xVel           = str(0.0)
nt             = str(10000)
outputEnabled  = str(1)
isSlipBoundary = str(0)
isPeriodic     = str(1)

### Demonstate:
* Energy stabilizes over time.
* Particle density stabilizes over time.
* Velocity profile after long time. 
'''


numSoParticles = str(2400)
Lx             = str(40)
Ly             = str(20)
Lz             = str(0)
Fe             = str(0.02)
xVel           = str(0.0)
nt             = str(40000)
outputEnabled  = str(1)
isSlipBoundary = str(0)
isPeriodic     = str(1)

for i in range(n): 
    subprocess.run(["./output/main", numSoParticles, Lx, Ly, Lz, Fe, xVel,
                   nt, outputEnabled, isSlipBoundary, isPeriodic, str(2020 + i)])
    if bool(outputEnabled): 
        os.rename('./outPos.csv', './data/outPosNVIDIATest4' +
                  str(i)+'.csv')
        os.rename('./outVel.csv', './data/outVelNVIDIATest4' +
                  str(i)+'.csv')