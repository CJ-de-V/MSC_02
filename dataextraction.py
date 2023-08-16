import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


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
basedir = os.getcwd()
df = pd.DataFrame(columns=["Np", "Nt", "Kp", "Kt", "lp", "Rg", "E2E", "zd", "lb", "sd_lp", "sd_Rg", "sd_zd"])
p1 = np.linspace(0., 0., 3)
p2 = np.linspace(0., 0., 3)
N = 0  # Number of particles
Nb = 0  # number of bonds
Nc = 0  # number of correlation functions
bondvectors = 0
avgbondlength = 0

Np = [8 , 32, 64, 128]  # [8 , 32, 64, 128] [0.25,  1, 4, 16, 64]
Nt = [8 , 32, 64, 128]
Kp = [0.25,  1, 4, 16, 64]
Kt = [0.25,  1, 4, 16, 64]

for nop in Np:
    for nt in Nt:
        for kp in Kp:
            for kt in Kt:
                discriminator = "Np_" + str(nop) + "_Nt_" + str(nt) + "_Kp_" + str(kp) + "_Kt_" + str(kt)
                os.chdir(basedir + "/output/" + discriminator)
                Nb = nop - 1
                Nc = Nb - 1
                bondvectors = np.zeros((Nb, 3))
                correlationfunctions = np.zeros(Nc)
                avgbondlength = 0
                e2edist = 0
                numavg = 0  # Number of timesteps we're averaging over, used for the normalization
                print('starting with analysis of: ' + discriminator)

                # persistence and bond length calculation
                with open('dump' + discriminator + '.dynamics') as datafile:
                    # Index-1 of the above indicates how many atoms separate the bondvectors
                    skiplines(datafile, 9)  # skips first boilerplate
                    line = 'liney'

                    while line != '':
                        line = (datafile.readline()).split(" ")
                        p1 = np.array([float(line[2]), float(line[3]), float(line[4])])
                        # position 1 1 is now real
                        p2 = np.array([0.0, 0.0, 0.0])
                        # reads particle 2 to N in and finds the bondvectors

                        startpos = p1 # position of leftmost bond
                        for i in range(1, nop):  # reading past entry 1 which we already read
                            line = (datafile.readline()).split(" ")
                            p2 = np.array([float(line[2]), float(line[3]), float(line[4])])
                            bondvectors[i - 1] = normalize(np.subtract(p2, p1))
                            p1 = p2
                            # if matches the bonding point record its Z position
                        E2Evector = np.subtract(p2,startpos)
                        e2edist = e2edist + np.sqrt(np.dot(E2Evector,E2Evector))

                        for i in range(Nc):  # iterates over different spacings particles can have
                            runninavg = 0.0
                            for j in range(0, Nb - i):  # iterates over all legal bonds with i bonds between them
                                runninavg += np.dot(bondvectors[j], bondvectors[j + i])  # Here be where we absed
                            correlationfunctions[i] += runninavg / (Nb - i)

                        skiplines(datafile, nt)  # skips the tether's lines
                        skiplines(datafile, 8)  # skips boilerplate
                        line = datafile.readline()  # reads in last line of boilerplate to confirm we're not at EOF
                        numavg += 1
                        # finished reading in and processing one timestep

                    avgbondlength = avgbondlength / (Nb * numavg)
                    e2edist=e2edist/numavg
                    y = correlationfunctions / numavg
                    x = np.arange(len(y)) * avgbondlength
                    # x range in LJ units
                    weights = np.reciprocal(np.arange(Nc + 0.0, 0.0, -1.0))
                    param, param_cov = curve_fit(test, x, y, maxfev=1000, sigma=weights)
                    lp = param[0]
                    sd_lp = np.sqrt(param_cov[0][0])

                # RG data
                with open('radius_of_gyration' + discriminator + '.dat') as datafile:
                    fulldat = datafile.readlines()  # these files are short enough for us to just use readlines
                    fulldat.pop(0)
                    for i in range(len(fulldat)):
                        fulldat[i] = float(fulldat[i].split(' ')[1])
                    Rg = np.mean(fulldat)
                    sd_Rg = np.var(fulldat)

                # Zdistance data

                with open('dump' + discriminator + '.dynamics') as datafile:
                    skiplines(datafile, 7)
                    zpeak = float(datafile.readline().split(" ")[1]) - 1
                    # finds walls Z pos, the -1 is to account for the fact that our wall is 1 unit away from the boundary

                with open('pcom' + discriminator + '.dat') as datafile:
                    fulldat = datafile.readlines()[2:]  # pops 2 first lines
                    for i in range(len(fulldat)):
                        fulldat[i] = float(fulldat[i].split(' ')[3])
                    numzdata = np.array(fulldat)
                    numzdata=zpeak-numzdata
                    zdis = np.mean(numzdata)
                    sd_zd = np.sqrt(np.var(numzdata))

                entry = [nop, nt, kp, kt, lp, Rg, e2edist , zdis,
                         avgbondlength, sd_lp, sd_Rg, sd_zd]
                datline = np.asarray(entry, dtype=float)
                df.loc[len(df)] = datline
print(df, '\n')
os.chdir(basedir)
df.to_csv('finaldata.csv')
# print you pandas here
