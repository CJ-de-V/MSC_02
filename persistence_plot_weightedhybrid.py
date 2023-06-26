from numpy.fft import fft, ifft
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import sys

# plots persistence length of the given dump file, graphical version of dataextraction.

avgbondlength = 0


def normalize(v):
    norm = np.linalg.norm(v)
    global avgbondlength
    avgbondlength += norm
    # this is only ever called to normalize bond vectors, but before that we simply
    # record the length of said vector, should possibly also average it sometime...
    if norm == 0:
        norm = np.finfo(v.dtype).eps
    return v / norm


def skiplines(f, n):
    for i in range(n):
        f.readline()

def test(x, b):
    return np.exp(-x / b)

np.set_printoptions(threshold=np.inf)
p1 = np.linspace(0., 0., 3)
p2 = np.linspace(0., 0., 3)
Ntot = 0 # total number of monomers
N = 0  # Number of monomers in POLYMER
Nt = 0 # Number of monomers in tether
Nb = 0  # number of bonds
Nc = 0  # number of correlation functions
nSkippableLines = 9
bondvectors = 0

with open(sys.argv[1]) as datafile:
     skiplines(datafile, 3)
     Ntot= int(datafile.readline())
     skiplines(datafile,5)
     for i in range(1,Ntot):
         if(datafile.readline().split(" ")[1]=="1"):
             N+=1

Nt=Ntot-N
Nb=N-1
Nc=Nb-1
bondvectors = np.zeros((Nb, 3))
correlationfunctions = np.zeros(Nc)

#printing all values here to confirm it's working properly
print("Ntot:" + str(Ntot)+ ", N:" + str(N)+ ", Ntet:" + str(Nt)+"\n")

line = 'liney'
numavg = 0  # Number of timesteps we're averaging over, used for the normalization
with open(sys.argv[1]) as datafile:
    skiplines(datafile, 9)  # skips boilerplate
    while line != '':
        # one run of this while loop innard is one timestep calculation
        line = (datafile.readline()).split(" ")
        print(line)
        p1 = np.array([float(line[2]), float(line[3]), float(line[4])])
        # position 1 1 is now real
        p2 = np.array([0.0, 0.0, 0.0])
        # reads particle 2-N in and finds the bondvectors
        for i in range(1, N):
            line = (datafile.readline()).split(" ")
            p2 = np.array([float(line[2]), float(line[3]), float(line[4])])
            bondvectors[i - 1] = normalize(np.subtract(p2, p1))
            p1 = p2

        # i is the separation between these vectors
        for i in range(Nc):  # iterates over different spacings particles can have
            runninavg = 0.0
            for j in range(0, Nb - i):  # iterates over all legal bonds with i bonds between them
                runninavg += np.dot(bondvectors[j], bondvectors[j + i])  # Here be where we absed
            correlationfunctions[i] += runninavg / (Nb - i)
        numavg += 1
        skiplines(datafile,Nt) # skips the tether's lines
        skiplines(datafile, 8) # skips boilerplate
        line=datafile.readline() #reads last line of boilerplate, if this is not present it means we're at EOF


avgbondlength = avgbondlength / (Nb * numavg)
y = correlationfunctions / numavg
x = np.arange(len(y)) * avgbondlength
# x range in LJ units
weights = np.reciprocal(np.arange(Nc+0.0, 0.0, -1.0))

# actual terminal dat
print('N: ' + str(N) + ', Nc:' + str(Nc) + ', Ntimesteps:' + str(numavg))

print(weights)
print(correlationfunctions)
# sigma is the relative errors in the y-data.... basically we have lower errors for low x since we have more data
param, param_cov = curve_fit(test, x, y, maxfev=1000, sigma=weights)
print("persistence length:")
print(param)
print("Covariance of coefficients:")
print(param_cov)
print("Average Bond Length: ")
print(avgbondlength)
plt.text(5, 0.9, 'PL: ' + str(param[0]) + 'with covar:' + str(param_cov[0, 0]) + '\n avg bond length: ' +
         str(avgbondlength), fontsize=14)
plt.xlim(left=0, right=max(x) + 1)
plt.plot(test(x, param[0]))
plt.plot(y, '.')

#plt.savefig(sys.argv[1].split('/')[1]+'.pdf')
plt.show()
