# batchcontrun.py
# Purpose: adapted from batchrun
# Should do the exact same thing as batchrun, except running from restart files is slightly different.

# Serial syntax: batchcont.py ???
#                ??????

# Parallel syntax: mpirun -np 4 python3 batchcont.py ???
#                  in.lammps = LAMMPS input script

from __future__ import print_function
import sys

import numpy as np
import pandas
import pandas as pd
import ctypes
from mpi4py import MPI
from lammps import lammps
import os
import pandas as pd

# basedirectory
basedir = os.getcwd()

me = 0
me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

N = [32, 64, 128, 256, 512]
K = [0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128]
Nruns = 500000

# List of N values along with which K values need to be run for a longer time.
# runmelonger = [[32, [0.25,0.5,1,2,4,8,16,32, 64, 128]],
#                [64, [16, 32, 64, 128]],
#                [128, [16]],
#                [256, [0.25, 0.5, 16]],
#                [512, [0.25, 0.5, 1, 16]]]


# adapt from continuedrestart?
# Sets up the computes and variables we require for our dump(s)
infile = 'in.restart'

for n in N:
    for k in K:
        os.chdir(basedir)  # back to file dir
        # change into dir where we wish to read/write... also using this since CD with multithreads to relative paths
        # is a bit wonky but this is functional
        os.chdir(basedir + '/N' + str(n))
        lmp = lammps()
        lmp.command('log log.N' + str(n) + 'k' + str(k))
        lmp.command('print \">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>N: ' + str(
            n) + ' K: ' + str(k) + '\"')

        # read in which restart file we're working with

        lmp.command('read_restart restart._N' + str(n) + '_k' + str(k) + '.dat')

        lines = open(basedir + '/' + infile, 'r').readlines()
        for line in lines:
            lmp.command(line)
        lmp.command("run 10000");
        # setup fixes that write to file, Thermo, R_g, and dump file, read in restart file
        fixsetup = ["fix mythermofile all print 10000 \"$t ${mytemp} ${myepair}\" file thermo_output_N" + str(
            n) + "_k" + str(k) + ".dat screen no",
                    # thermodynamic data outputted to the file appropriately named (sort of)
                    "fix myRG2file all print 10000 \"$t ${RG2}\" file radius_of_gyration_squared_N" + str(
                        n) + "_k" + str(k) + ".dat screen no",
                    "dump dum2 all custom 10000 dump_N" + str(n) + "_k" + str(k) + ".dynamics id type x y z",
                    "dump_modify dum2  sort id"
                    # ,"fix recenterinator all recenter 0.0 0.0 0.0"
                    ]
        lmp.commands_list(fixsetup)

        # run and collect data under these new fixes
        lmp.command("run " + str(Nruns))

        # End of SimUlation Events
        cleansetup = ["write_restart 	restart._N" + str(n) + "_k" + str(k) + ".dat"]
        lmp.commands_list(cleansetup)
        lmp.close()
# uncomment if running in parallel via mpi4py
print("Proc %d out of %d procs has" % (me, nprocs), lmp)

MPI.COMM_WORLD.Barrier()
# Definitely worth looking at the for throwing these in the for loops, might be that only one core is ending the
# dumps before the others get to dumping
MPI.Finalize()
