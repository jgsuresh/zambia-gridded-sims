"""
Run multi-node simulation
"""

# Calls setup_run.py

# import os
import subprocess

# Run mode 0 builds input files for simulations
# Run mode 1 submits simulations based on these input files.
run_mode = 0

# milen_catch_list = ['chisanga', 'sinamalima', 'chiyabi', 'sinafala', 'munyumbwe', 'chipepo',
#                     'chabbobboma', 'bbondo', 'luumbo', 'nyanga_chaamwe']
milen_catch_list = ['bbondo','chabbobboma','chisanga','chiyabi','luumbo','munyumbwe','nyanga chaamwe','sinafala','sinamalima']
priority = "Highest"
coreset = "emod_abcd"
num_cores = 12

#######################################################################################

arg_dict = {}
arg_dict["Normal"] = 0
arg_dict["AboveNormal"] = 1
arg_dict["Highest"] = 2
arg_dict["emod_32cores"] = 0
arg_dict["emod_abcd"] = 1


for catch in milen_catch_list:
    cmd = "python setup_run.py {} {} {} {} {}".format(run_mode,arg_dict[priority],arg_dict[coreset],num_cores,catch)
    print("Running {}".format(cmd))
    # os.system(cmd)
    # subprocess.call(["python","setup_run.py",run_mode,arg_dict[priority],arg_dict[coreset],num_cores,catch])
    subprocess.call(cmd)



# Process input:
# argv[1] sets whether the script will create simulation input files, or submit simulation jobs to COMPS

# argv[2] sets run priority
# if sys.argv[2] == 0:
#     priority = "Normal"
# elif sys.argv[2] == 1:
#     priority = "AboveNormal"
# elif sys.argv[2] == 2:
#     priority = "Highest"
#
# # argv[3] sets which coreset the job will run on.  emod_32cores is faster but should only be used for test runs
# if sys.argv[3] == 0:
#     coreset = "emod_abcd"
# elif sys.argv[3] == 1:
#     coreset = "emod_32cores"
#
# # argv[4] sets number of cores to run job on
# num_cores = sys.argv[4]
#
# # arg[5] sets the catchment name:
# catch = sys.argv[5]