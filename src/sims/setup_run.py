"""
Run multi-node simulation
"""

from experiment_setup import *
from simtools.SetupParser import SetupParser
import sys

# ===================================================================================
try:
    # Process input:
    # argv[1] sets whether the script will create simulation input files, or submit simulation jobs to COMPS
    run_mode = sys.argv[1]

    # argv[2] sets run priority
    if sys.argv[2] == 0:
        priority = "Normal"
    elif sys.argv[2] == 1:
        priority = "AboveNormal"
    elif sys.argv[2] == 2:
        priority = "Highest"

    # argv[3] sets which coreset the job will run on.  emod_32cores is faster but should only be used for test runs
    if sys.argv[3] == 0:
        coreset = "emod_abcd"
    elif sys.argv[3] == 1:
        coreset = "emod_32cores"

    # argv[4] sets number of cores to run job on
    num_cores = sys.argv[4]

    # arg[5] sets the catchment name:
    catch = sys.argv[5]

except:
    print("Incorrect or missing inputs ")
    run_mode = 1
    priority = "Normal"
    coreset = "emod_abcd"
    num_cores = 12
    milen_catch_list = ['bbondo', 'chabbobboma', 'chisanga', 'chiyabi', 'luumbo', 'munyumbwe', 'nyanga chaamwe',
                        'sinafala', 'sinamalima']
    catch = milen_catch_list[6]


# ===================================================================================

if '_' in catch:
    catch = catch.replace('_',' ')
exp_name = '{}_mod_hab_params_v1'.format(catch)

gravity_migr_params = np.array([7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00])

base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
grid_pop_csv_file = base + 'data/gridded_pop/cleaned/all_max_pop.csv' # Max pop in grid cell for entire region
imm_1node_fp = base + "data/immunity/Sinamalima_1_node_immune_init_p1_33_p2_117.json"

# Intervention file names:
healthseek_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_healthseek_events.csv'
itn_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_itn_events.csv'
irs_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_irs_events.csv'
msat_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_msat_events.csv'
mda_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_mda_events.csv'
stepd_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_stepd_events.csv'


start_year = 2007
num_years = 12

comps_exp = COMPS_Experiment(base,
                             exp_name,
                             catch=catch,
                             grid_pop_csv_file=grid_pop_csv_file,
                             imm_1node_fp=imm_1node_fp,
                             migration_on=True,
                             start_year=start_year,
                             sim_length_years=num_years,
                             rcd_people_num=10,
                             gravity_migr_params=gravity_migr_params,
                             num_cores=num_cores,
                             healthseek_fn=healthseek_fn,
                             itn_fn=itn_fn,
                             irs_fn=irs_fn,
                             msat_fn=msat_fn,
                             mda_fn=mda_fn,
                             stepd_fn=stepd_fn,
                             larval_params_mode="milen",
                             immunity_mode="milen",
                             fudge_milen_habitats=True)

if run_mode == 0:
    comps_exp.file_setup(generate_climate_files=False)  # 1. Run this first to set up input files.  When submitting experiment to COMPS, comment this line out.

if __name__ == "__main__":

    if run_mode == 1:
        # 2. Run these files after input files are set up.  This will not work correctly if the file_setup() line is not commented out, sorry!
        SetupParser.init()

        SetupParser.set("HPC", "priority", priority)
        SetupParser.set("HPC", "node_group", coreset)


        comps_exp.submit_experiment(num_seeds=4,simple_intervention_sweep=False)
