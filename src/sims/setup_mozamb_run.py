"""
Run multi-node simulation
"""

from mozambique_experiments import MozambiqueExperiment
from experiment_setup import GriddedInputFilesCreator
from simtools.SetupParser import SetupParser


run_mode = 1
priority = "BelowNormal"
coreset = "emod_abcd"
parser_location = "HPC"
num_cores = 12
mozamb_catch_list = ["Caputine","Mahel","Panjane","Mapulanguene","Chicutso","Moine","Motaze","Chichuco","Facazissa","Magude-Sede"]
catch = mozamb_catch_list[0]


# ===================================================================================

exp_name = '{}_selfc_test'.format(catch)
start_year = 2007
num_years = 12

base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
grid_pop_csv_file = base + 'data/mozambique/grid_population.csv' # Max pop in grid cell for entire region

# Intervention file names:
healthseek_fn = base + 'data/mozambique/grid_all_healthseek_events.csv'
itn_fn = base + 'data/mozambique/grid_all_itn_events.csv'
irs_fn = base + 'data/mozambique/grid_all_irs_events.csv'
msat_fn = None
mda_fn = base + 'data/mozambique/grid_all_mda_events.csv'
stepd_fn = None


# Build config-builder:
mozamb_exp = MozambiqueExperiment(base,
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
                                            mozamb_exp.desired_cells,
                                            mozamb_exp.cb,
                                            grid_pop_csv_file,
                                            region=mozamb_exp.region,
                                            start_year=start_year,
                                            sim_length_years=num_years,
                                            immunity_mode="naive",
                                            larval_param_func=mozamb_exp.larval_params_func_for_calibration
                                            )

if __name__ == "__main__":

    if run_mode == 1:
        # 2. Run these files after input files are set up.  This will not work correctly if the file_setup() line is not commented out, sorry!
        SetupParser.init()

        SetupParser.set("HPC", "priority", priority)
        SetupParser.set("HPC", "node_group", coreset)

        mozamb_exp.submit_experiment(num_seeds=4)
