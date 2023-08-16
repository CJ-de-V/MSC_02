# Batchrun.py: runs a batch of LAMMPS experiments for the tethered version
# arguments: NONE
# requirements: in.setup, in.continue and a.out should be present in the run directory

import os
from mpi4py import MPI
from lammps import lammps
from subprocess import call

me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
basedir = os.getcwd()
Np = [8 , 32, 64, 128 , 256]  # [8 ,16, 32, 64, 128] [0.25,  1, 4, 16, 64]
Nt = [8 , 32, 64,128 , 256]
Kp = [0.25,  1, 4, 16, 64]
Kt = [0.25,  1, 4, 16, 64]
restarting = True

for np in Np:
    for nt in Nt:
        for kp in Kp:
            for kt in Kt:
                discriminator = "Np_" + str(np) + "_Nt_" + str(nt) + "_Kp_" + str(kp) + "_Kt_" + str(
                    kt)
                workingdir = basedir + "/output/" + discriminator
                os.makedirs(workingdir, exist_ok=True)
                os.chdir(workingdir)

                if not os.path.isfile(workingdir + "/config.tether"):
                    # print("no starting configuration, making one")
                    call(["../../a.out", str(np), str(nt), "config.tether"])

                if not os.path.isfile(workingdir + "/restart." + discriminator):
                    # print("no restart file, going to run sim from scratch")
                    # No restart file
                    restarting = False

                lmp = lammps()

                restartorfresh = ""
                if restarting:
                    lmp.command('read_restart restart.' + discriminator)
                    lines = open(basedir + '/' + "in.restart", 'r').readlines()
                    for line in lines:
                        lmp.command(line)
                    restartorfresh=restartorfresh+"Restart"

                else:
                    lines = open(basedir + '/' + "in.setup", 'r').readlines()
                    for line in lines:
                        lmp.command(line)
                    restartorfresh = restartorfresh + "Fresh"

                rundetails = "print " + "\"" + ">>>>>>>>>>>>>>>>DETAILS OF RUN<<<<<<<<<<<<<<<<" +"\n>>>>>>>>>>>>>>>>" + discriminator +"___"+restartorfresh+ "<<<<<<<<<<<<<<<<\""
                lmp.command(rundetails)
                # Commands needed for both setup and restart runs
                primersetup = ['log log.' + discriminator,
                               "fix mythermofile all print 10000 \"$t ${mytemp} ${myepair}\" file thermo_output" + discriminator + ".dat screen no",
                               # thermodynamic data outputted to the file appropriately named (sort of)
                               "fix myRG2file all print 10000 \"$t ${RG2}\" file radius_of_gyration" + discriminator + ".dat screen no",
                               "fix comfix polymer ave/time 1 1 10000 c_mycomcompute[*] file pcom"+discriminator+".dat",
                               #every 10 000 timesteps spits out COM of polymer
                               "dump dum2 all custom 10000 dump" + discriminator + ".dynamics id type x y z",
                               "dump_modify dum2  sort id",
                               "angle_coeff   1  " + str(kp),
                               "angle_coeff   2  " + str(kt)
                               ]
                lmp.commands_list(primersetup)

                lmp.command("run 100000")

                cleansetup = ["write_restart restart." + discriminator]
                lmp.commands_list(cleansetup)
                lmp.close()

print("Proc %d out of %d procs has" % (me, nprocs), lmp)
