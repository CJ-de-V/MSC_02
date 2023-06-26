import numpy
import numpy as np
from scipy.optimize import curve_fit

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


datfile = open("data.csv", "w")

np.set_printoptions(threshold=np.inf)
p1 = np.linspace(0., 0., 3)
p2 = np.linspace(0., 0., 3)
N = 0  # Number of particles
Nb = 0  # number of bonds
Nc = 0  # number of correlation functions
bondvectors = 0

Nmonomers = [32, 64, 128, 256, 512]
kvalues = [0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128]
datfile.write("N, k, lp, lp_sd, R_g, R_g_sd, lb\n")
for nmonom in Nmonomers:
    print('starting with N=' + str(nmonom))
    for kval in kvalues:
        with open('N' + str(nmonom) + '/dump_N' + str(nmonom) + '_k' + str(kval) + '.dynamics') as datafile:
            avgbondlength = 0
            skiplines(datafile, 3)
            N = int(datafile.readline())
            Nb = N - 1
            skiplines(datafile, 5)
            bondvectors = np.zeros((Nb, 3))

            Nc = Nb - 1
            correlationfunctions = np.zeros(Nc)
            # Index-1 of the above indicates how many atoms separate the bondvectors

            line = 'liney'
            numavg = 0  # Number of timesteps we're averaging over, used for the normalization
            while line != '':
                line = (datafile.readline()).split(" ")
                p1 = np.array([float(line[2]), float(line[3]), float(line[4])])
                # position 1 1 is now real
                p2 = np.array([0.0, 0.0, 0.0])
                # reads particle 2-N in and finds the bondvectors
                for i in range(1, N):
                    line = (datafile.readline()).split(" ")
                    p2 = np.array([float(line[2]), float(line[3]), float(line[4])])
                    bondvectors[i - 1] = normalize(np.subtract(p2, p1))
                    p1 = p2

                for i in range(Nc):  # iterates over different spacings particles can have
                    runninavg = 0.0
                    for j in range(0, Nb - i):  # iterates over all legal bonds with i bonds between them
                        runninavg += np.dot(bondvectors[j], bondvectors[j + i])  # Here be where we absed
                    correlationfunctions[i] += runninavg / (Nb - i)
                skiplines(datafile, 8)
                line = datafile.readline()
                numavg += 1
                # finished reading in and processing one timestep

            avgbondlength = avgbondlength / (Nb * numavg)
            y = correlationfunctions / numavg
            x = np.arange(len(y)) * avgbondlength
            # x range in LJ units
            weights = np.reciprocal(np.arange(Nc + 0.0, 0.0, -1.0))
            param, param_cov = curve_fit(test, x, y, maxfev=1000, sigma=weights)
            lp = param[0]
            # R_gs hootenanny

            with open('N' + str(nmonom) + '/radius_of_gyration_squared_N' + str(nmonom) + '_k' + str(
                    kval) + '.dat') as datafile:
                # readandmakebusiness
                fulldat = datafile.readlines()
                fulldat.pop(0)
                for i in range(len(fulldat)):
                    fulldat[i] = float(fulldat[i].split(' ')[1])
                Rg = numpy.mean(fulldat)
                RgVar = numpy.var(fulldat)

            datfile.write(str(nmonom) + ', ')  # N
            datfile.write(str(kval) + ', ')  # k
            datfile.write(str(lp) + ', ')  # lp
            # below is from switching to LJ units from bondvector units
            datfile.write(str(np.sqrt(param_cov[0][0])) + ', ')  # SD of lp

            datfile.write(str(Rg) + ', ')  # Rg
            datfile.write(str(np.sqrt(RgVar)) + ', ')  # SD of Rg
            datfile.write(str(avgbondlength) + '\n')  # lb
# datfile.close()
