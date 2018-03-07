# Goal: do a simple 2D calibration where later runs draw from serialized files of best-fit earlier round.
# Round 0: Long burn-in, and serialize files before intervention (requires serialize_year_first_round)
# Round N>0: (need a start year) use serialized files of best fit from previous run.  Start it 5 years before the serialized_year_first_round



"""
Run multi-node simulation
"""

import pandas as pd
import copy
import os.path

from calibtool.CalibManager import CalibManager
from calibtool.algorithms.OptimTool import OptimTool
from mozambique_experiments import MozambiqueExperiment
from experiment_setup import GriddedInputFilesCreator
from simtools.SetupParser import SetupParser
from dtk.interventions.habitat_scale import scale_larval_habitats

from GriddedCalibSite import GriddedCalibSite
from calibtool.plotters.LikelihoodPlotter import LikelihoodPlotter
from calibtool.plotters.OptimToolPlotter import OptimToolPlotter
from calibtool.plotters.SiteDataPlotter import SiteDataPlotter
from dtk.vector.species import set_larval_habitat

run_mode = 1
priority = "Normal"
coreset = "emod_abcd"
parser_location = "HPC"
num_cores = 8
mozamb_catch_list = ["Caputine","Mahel","Panjane","Mapulanguene","Chicutso","Moine","Motaze","Chichuco","Facazissa","Magude-Sede"]
catch = mozamb_catch_list[0]

COMPS_calib_exp_name = 'Serial_{}'.format(catch)

bairro_level_calibration = False
# ===================================================================================

exp_name = '{}_serial_test'.format(catch)

end_year = 2020
serialize_year = 2010
nonburnin_sim_length_years = 11
burnin_sim_length_years = 11 #60
nonburnin_sim_start_year = end_year - nonburnin_sim_length_years
burnin_sim_start_year = end_year - burnin_sim_length_years

iter0_write_time = 365 * (serialize_year - burnin_sim_start_year)
iterN_write_time = 365 * (serialize_year - nonburnin_sim_start_year)


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
                              start_year=burnin_sim_start_year,
                              sim_length_years=burnin_sim_length_years,
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
                                            start_year=burnin_sim_start_year,
                                            sim_length_years=burnin_sim_length_years,
                                            immunity_mode="naive",
                                            larval_param_func=mozamb_exp.larval_params_func_for_calibration
                                            )





# Calibration-specific stuff:
sites = [GriddedCalibSite()]
# The default plotters used in an Optimization with OptimTool
plotters = [LikelihoodPlotter(combine_sites=True),
            SiteDataPlotter(num_to_plot=5, combine_sites=True),
            OptimToolPlotter()  # OTP must be last because it calls gc.collect()
            ]

# params = [
#     {
#         'Name': 'funestus_scale',
#         'Dynamic': True,
#         'MapTo': 'funestus_scale', # <-- DEMO: Custom mapping, see map_sample_to_model_input below
#         'Guess': 8.5,
#         'Min': 6,
#         'Max': 9.5
#     },
#     {
#         'Name': 'arabiensis_scale',
#         'Dynamic': True,
#         'MapTo': 'arabiensis_scale',
#         'Guess': 9,
#         'Min': 6,
#         'Max': 9.5
#     },
# ]

if not bairro_level_calibration:
    params = [
        {
            'Name': 'arabiensis_scale',
            'Dynamic': True,
            'MapTo': 'arabiensis_scale',
            'Guess': 8.7,
            'Min': 8.2,
            'Max': 9.1
        },
        {
            'Name': 'arabiensis_funestus_ratio',
            'Dynamic': True,
            'MapTo': 'arabiensis_funestus_ratio',  # <-- DEMO: Custom mapping, see map_sample_to_model_input below
            'Guess': 5.,
            'Min': 3.,
            'Max': 12.
        },
    ]
else:
    # How many bairros?  2x params for each.
    bairro_dict = MozambiqueExperiment.find_bairros_for_this_catchment(catch)
    n_bairros = bairro_dict["num_bairros"]

    # Based on # of params, set up parameter dictionary
    params = [
        {
            'Name': 'arabiensis_scale',
            'Dynamic': True,
            'MapTo': 'arabiensis_scale',
            'Guess': 8.7,
            'Min': 8.2,
            'Max': 9.1
        },
        {
            'Name': 'arabiensis_funestus_ratio',
            'Dynamic': True,
            'MapTo': 'arabiensis_funestus_ratio',  # <-- DEMO: Custom mapping, see map_sample_to_model_input below
            'Guess': 5.,
            'Min': 3.,
            'Max': 12.
        }
    ]

    for bairro_num in bairro_dict["bairo_num_list"]:
        grid_cells_this_bairro = bairro_dict[bairro_num]
        params_this_barrio = [
            {
                'Name': 'b{}_arab_mult'.format(bairro_num),
                'Dynamic': True,
                'MapTo': 'b{}_arab_mult'.format(bairro_num),
                'Guess': 1.0,
                'Min': 0.333,
                'Max': 3.
            },
            {
                'Name': 'b{}_funest_mult'.format(bairro_num),
                'Dynamic': True,
                'MapTo': 'b{}_funest_mult'.format(bairro_num),
                'Guess': 1.0,
                'Min': 0.333,
                'Max': 3.
            }
        ]

        params.append(params_this_barrio)



def map_sample_to_model_input(cb, sample):
    # Serialization:
    # Check if can find LL_all.csv:
    LL_all_path = COMPS_calib_exp_name + "/_plots/LL_all.csv"


    # If you can't find LL_all, then we are in the zeroth iteration
    if not os.path.isfile(LL_all_path):
        # Do a long burn-in run, and serialize files
        # serialization_read_time = 50 * 365
        start_year = burnin_sim_start_year
        sim_length_years = burnin_sim_length_years
        serialization_write_time = iter0_write_time


    # If you can find LL_all, then we are not in the zeroth iteration right now.
    else:
        # If yes, then use it to find best fit so far.  Get output path for that run.  Set the serialization path to be there
        LL_all = pd.read_csv(LL_all_path)
        dir_list = list(LL_all['outputs'])
        best_run_dir = dir_list[-1]
        # If >1 seeds per parameter set, then this will be a comma separated list.  Just take the first in the list
        # (Thankfully this line will run even if there is no comma in the string)
        best_run_dir = best_run_dir.split(',')[0]

        start_year = nonburnin_sim_start_year
        sim_length_years = nonburnin_sim_length_years
        # Find out which iteration we are on.
        # change start year, and run length
        best_iteration = list(LL_all["iteration"])[-1]
        if best_iteration == 0:
            serialization_read_time = iter0_write_time
            serialization_write_time = iterN_write_time
        else:
            serialization_read_time = iterN_write_time
            serialization_write_time = iterN_write_time


    mozamb_exp.start_year = start_year
    mozamb_exp.sim_length_years = sim_length_years
    mozamb_exp.implement_baseline_healthseeking(cb)
    mozamb_exp.implement_interventions(cb,True,True,False,True,False)

    # Add one crazy year with 3x larval habitats
    scale_larval_habitats(cb,
                          pd.DataFrame({'LINEAR_SPLINE': [3]}),
                          start_day=365 * (2017 - start_year))
    scale_larval_habitats(cb,
                          pd.DataFrame({'LINEAR_SPLINE': [1]}),
                          start_day=365 * (2018 - start_year))


    # Now that cb has interventions added, give it the needed serialization information:
    cb.update_params({"Serialization_Time_Steps": [serialization_write_time]})
    if os.path.isfile(LL_all_path):
        cb.update_params({
            "Serialized_Population_Path": best_run_dir + "/output",
            'Serialized_Population_Filenames': ['state-%05d-%03d.dtk' % (serialization_read_time, corenum) for corenum in
                                                range(num_cores)]
            # 'Serialized_Population_Filenames': ['state-%03d.dtk' % corenum for corenum in range(num_cores)]
        })



    a_sc = sample['arabiensis_scale']
    arab_funest_ratio = sample['arabiensis_funestus_ratio']

    # Prashanth ento
    hab = {
        'arabiensis': {
            "TEMPORARY_RAINFALL": 2e3*pow(10,a_sc-8),
            "LINEAR_SPLINE": {
                "Capacity_Distribution_Per_Year": {
                    "Times": [0.0, 30.417, 60.833, 91.25, 121.667, 152.083, 182.5, 212.917, 243.333, 273.75, 304.167,
                              334.583],
                    "Values": [0.273953355,
                               4.226011848,
                               5.140191814,
                               9.363408701,
                               0.0,
                               0.414082115,
                               0.139915067,
                               0.186456901,
                               0.015611024,
                               0.101027567,
                               0.0,
                               0.121014426
                               ]
                },
                "Max_Larval_Capacity": pow(10,a_sc)
            }
        },
        'funestus': {
            "WATER_VEGETATION": 2.2e7*pow(10,a_sc-8)/arab_funest_ratio,
            "LINEAR_SPLINE": {
                "Capacity_Distribution_Per_Year": {
                    "Times": [0.0, 30.417, 60.833, 91.25, 121.667, 152.083, 182.5, 212.917, 243.333, 273.75, 304.167,
                              334.583],
                    # "Values": [0.0, 0.0, 0.0, 0.2, 0.8, 1.0, 1.0, 1.0, 0.5, 0.2, 0.0, 0.0] # Caitlin
                    # "Values": [0.0, 0.0, 0.0, 0.0, 0.2, 1.0, 1.0, 1.0, 0.5, 0.2, 0.0, 0.0] # Milen
                    "Values": [0.0,
                               1.202730029,
                               0.112447779,
                               1.467850365,
                               2.470962168,
                               1.064668156,
                               4.806719314,
                               0.914212162,
                               9.919572963,
                               0.437353893,
                               0.392657387,
                               1.213697659
                               ]
                },
                "Max_Larval_Capacity": pow(10,a_sc)/arab_funest_ratio
            }
        }
    }

    set_larval_habitat(cb, hab)

    tags = sample
    return tags


optimtool = OptimTool(params,
                      samples_per_iteration=18,
                      center_repeats=1,
                      sigma_r=0.1) #increase radius of hypersphere (default is sigma_r=0.02)

foo = copy.deepcopy(mozamb_exp.cb_no_interventions)

calib_manager = CalibManager(name=COMPS_calib_exp_name,
                             config_builder=mozamb_exp.cb,
                             map_sample_to_model_input_fn=map_sample_to_model_input,
                             sites=sites,
                             next_point=optimtool,
                             sim_runs_per_param_set=5,
                             max_iterations=50,
                             plotters=plotters)

run_calib_args = {
    "calib_manager": calib_manager
}


if __name__ == "__main__":

    if run_mode == 1:
        # 2. Run these files after input files are set up.  This will not work correctly if the file_setup() line is not commented out, sorry!

        if parser_location == "LOCAL":
            SetupParser.init("LOCAL")

        else:
            SetupParser.init()

            SetupParser.set("HPC", "priority", priority)
            SetupParser.set("HPC", "node_group", coreset)

        cm = run_calib_args["calib_manager"]
        cm.run_calibration()

        # comps_exp.submit_experiment(num_seeds=4,simple_intervention_sweep=False)
