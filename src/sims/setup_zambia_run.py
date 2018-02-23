"""
Run multi-node simulation
"""

from zambia_experiments import ZambiaExperiment
from experiment_setup import GriddedInputFilesCreator
from simtools.SetupParser import SetupParser


run_mode = 1
priority = "BelowNormal"
coreset = "emod_abcd"
parser_location = "HPC"
num_cores = 12
milen_catch_list = ['bbondo', 'chabbobboma', 'chisanga', 'chiyabi', 'luumbo', 'munyumbwe', 'nyanga chaamwe',
                    'sinafala', 'sinamalima']
catch = milen_catch_list[0]


# ===================================================================================

exp_name = '{}_selfc_test'.format(catch)
start_year = 2007
num_years = 12

base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
grid_pop_csv_file = base + 'data/zambia/cleaned/all_max_pop.csv' # Max pop in grid cell for entire region

# Intervention file names:
healthseek_fn = base + 'data/zambia/grid_all_healthseek_events.csv'
itn_fn = base + 'data/zambia/grid_all_itn_events.csv'
irs_fn = base + 'data/zambia/grid_all_irs_events.csv'
msat_fn = base + 'data/zambia/grid_all_msat_events.csv'
mda_fn = base + 'data/zambia/grid_all_mda_events.csv'
stepd_fn = base + 'data/zambia/grid_all_stepd_events.csv'


# Build config-builder:
zambia_exp = ZambiaExperiment(base,
                              exp_name,
                              catch,
                              healthseek_fn=healthseek_fn,
                              itn_fn=itn_fn,
                              irs_fn=irs_fn,
                              msat_fn=msat_fn,
                              mda_fn=mda_fn,
                              stepd_fn=stepd_fn,
                              start_year=start_year,
                              sim_length_years=num_years,
                              immunity_mode="naive",
                              num_cores=num_cores,
                              parser_location=parser_location)

# Create necessary input files
if run_mode == 0:
    file_creator = GriddedInputFilesCreator(base,
                                            exp_name,
                                            zambia_exp.desired_cells,
                                            zambia_exp.cb,
                                            grid_pop_csv_file,
                                            region=zambia_exp.region,
                                            start_year=start_year,
                                            sim_length_years=num_years,
                                            immunity_mode="naive",
                                            larval_param_func=zambia_exp.larval_params_func_milen
                                            )

if __name__ == "__main__":

    if run_mode == 1:
        # 2. Run these files after input files are set up.  This will not work correctly if the file_setup() line is not commented out, sorry!
        SetupParser.init()

        SetupParser.set("HPC", "priority", priority)
        SetupParser.set("HPC", "node_group", coreset)


        zambia_exp.submit_experiment(num_seeds=4)
