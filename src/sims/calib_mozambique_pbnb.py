"""
Calibrate multi-node simulation
"""

import copy

from experiment_setup import *
from simtools.SetupParser import SetupParser

from calibtool.CalibManager import CalibManager
from calibtool.algorithms.OptimTool import OptimTool
from calibtool.plotters.LikelihoodPlotter import LikelihoodPlotter
from calibtool.plotters.OptimToolPlotter import OptimToolPlotter
from calibtool.plotters.SiteDataPlotter import SiteDataPlotter
from simtools.AssetManager.SimulationAssets import SimulationAssets
from dtk.vector.species import set_larval_habitat
from dtk.interventions.habitat_scale import scale_larval_habitats

from calibtool.algorithms.PBnB.OptimTool_PBnB import OptimTool_PBnB, par
# from calibtool.algorithms.PBnB.OptimTool_PBnB import OptimTool_PBnB, par
from calibtool.plotters.OptimToolPBnBPlotter import OptimToolPBnBPlotter

from dtk.vector.study_sites import configure_site

from GriddedCalibSite import GriddedCalibSite

import matplotlib
matplotlib.use('Agg')


run_mode = 1
run_location = "HPC" #""LOCAL"
priority = "AboveNormal"
coreset = "emod_abcd" #"emod_32cores"
num_cores = 8

mozambique_catch_list = ["Caputine","Mahel","Panjane","Mapulanguene","Chicutso","Moine","Motaze","Chichuco","Facazissa","Magude-Sede"]
catch = mozambique_catch_list[0]

# ===================================================================================

exp_name = '{}_calibtest'.format(catch)

gravity_migr_params = np.array([7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00])

base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
grid_pop_csv_file = base + 'data/mozambique/grid_population.csv' # Max pop in grid cell for entire region

# Intervention file names:
healthseek_fn = base + 'data/mozambique/grid_all_healthseek_events.csv'
itn_fn = base + 'data/mozambique/grid_all_itn_events.csv'
irs_fn = base + 'data/mozambique/grid_all_irs_events.csv'
# msat_fn = base + 'data/mozambique/grid_all_msat_events.csv'
mda_fn = base + 'data/mozambique/grid_all_mda_events.csv'
# stepd_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_stepd_events.csv'


# start_year = 2000
# num_years = 17
# start_year = 2010
# num_years = 0.03 #1
start_year = 1994
num_years = 24


comps_exp = COMPS_Experiment(base,
                             exp_name,
                             catch=catch,
                             grid_pop_csv_file=grid_pop_csv_file,
                             migration_on=True,
                             start_year=start_year,
                             sim_length_years=num_years,
                             # rcd_people_num=10,
                             gravity_migr_params=gravity_migr_params,
                             num_cores=num_cores,
                             healthseek_fn=healthseek_fn,
                             itn_fn=itn_fn,
                             irs_fn=irs_fn,
                             # msat_fn=msat_fn,
                             mda_fn=mda_fn,
                             # stepd_fn=stepd_fn,
                             larval_params_mode="calibrate",
                             immunity_mode="naive")

if run_mode == 0:
    comps_exp.file_setup()  # 1. Run this first to set up input files.  When submitting experiment to COMPS, comment this line out.

cb = comps_exp.return_cb_for_calibration()


# Prashanth:
cb.update_params({'Antigen_Switch_Rate': pow(10, -9.116590124),
                  'Base_Gametocyte_Production_Rate': 0.06150582,
                  'Base_Gametocyte_Mosquito_Survival_Rate': 0.002011099,

                  'Falciparum_MSP_Variants': 32,
                  'Falciparum_Nonspecific_Types': 76,
                  'Falciparum_PfEMP1_Variants': 1070,
                  'Gametocyte_Stage_Survival_Rate': 0.588569307,

                  'MSP1_Merozoite_Kill_Fraction': 0.511735322,
                  'Max_Individual_Infections': 3,
                  'Nonspecific_Antigenicity_Factor': 0.415111634,
                  })


# Add specific boosts to larval parameters

# df = pd.DataFrame( { 'NodeID' : [1001, 1, 2],
#                      'CONSTANT': [0, 0, 1],
#                      'TEMPORARY_RAINFALL': [0, 1, 0],
#                      'species': ['arabiensis', 'arabiensis', 'funestus']
#                      })


# Add one crazy year with 3x larval habitats
scale_larval_habitats(cb,
                      pd.DataFrame({'LINEAR_SPLINE': [3]}),
                      start_day=365*(2017-start_year))
scale_larval_habitats(cb,
                      pd.DataFrame({'LINEAR_SPLINE': [1]}),
                      start_day=365*(2018-start_year))
# scale_larval_habitats(cb,
#                       pd.DataFrame({'LINEAR_SPLINE': [3],
#                                     'Start_Day': [365*(2018-start_year)]}),
#                       start_day=365*(2018-start_year))


# Calibration-specific stuff:
sites = [GriddedCalibSite()]
# The default plotters used in an Optimization with OptimTool
# plotters = [LikelihoodPlotter(combine_sites=True),
#             SiteDataPlotter(num_to_plot=5, combine_sites=True),
#             OptimToolPlotter()  # OTP must be last because it calls gc.collect()
#             ]
plotters = [LikelihoodPlotter(combine_sites=True),
            SiteDataPlotter(num_to_plot=5, combine_sites=True),
            OptimToolPBnBPlotter()]

params = [
    {
        'Name': 'funestus_scale',
        'Dynamic': True,
        'MapTo': 'funestus_scale', # <-- DEMO: Custom mapping, see map_sample_to_model_input below
        # 'Guess': 25,
        # 'Min': 0.1,
        # 'Max': 150
        # 'Guess': 6,
        # 'Min': 4,
        # 'Max': 8
        'Guess': 8.5,
        'Min': 7,
        'Max': 9.5#10
    },
    {
        'Name': 'arabiensis_scale',
        'Dynamic': True,
        'MapTo': 'arabiensis_scale',
        'Guess': 8.5,
        'Min': 6.5,
        'Max': 10.5
    },
]


def map_sample_to_model_input(cb, sample):
    tags = {}
    sample = copy.deepcopy(sample)

    if 'arabiensis_scale' in sample:
        a_sc = sample.pop('arabiensis_scale')

        # Zambia
        # hab = {'arabiensis': {'TEMPORARY_RAINFALL': 1e8 * a_sc, 'CONSTANT': 2e6}}

        # Prashanth ento
        hab = {'arabiensis': {
            # "TEMPORARY_RAINFALL": 2.2e7,
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
                # "Max_Larval_Capacity": 1e8 * f_sc
                "Max_Larval_Capacity": pow(10,a_sc)
            }
        }
        }

        set_larval_habitat(cb, hab)

        tags.update({'arabiensis_scale': a_sc})

    if 'funestus_scale' in sample:
        f_sc = sample.pop('funestus_scale')

        hab = {'funestus': {
            # "WATER_VEGETATION": 2e3, # Prashanth    #2e6, # Caitlin.  Milen had 2e7,
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
                "Max_Larval_Capacity": pow(10,f_sc)  # Prashanth
                # "Max_Larval_Capacity": 1e8 * f_sc # older
            }
        }
        }

        set_larval_habitat(cb, hab)

        tags.update({'funestus_scale': f_sc})

    for name,value in sample.items():
        print('UNUSED PARAMETER:', name)
    assert( len(sample) == 0 ) # All params used

    return tags









name = 'TestPBNB_{}'.format(catch)
#
# optimtool = OptimTool(params,
#                       samples_per_iteration=25,
#                       center_repeats=1)
optimtool_PBnB = OptimTool_PBnB(params,
                                s_running_file_name=name,
                                s_problem_type="deterministic",  # deterministic or noise
                                f_delta=par.f_delta,  # <-- to determine the quantile for the target level set
                                f_alpha=par.f_alpha,  # <-- to determine the quality of the level set approximation
                                i_k_b=par.i_k_b,  # <-- maximum number of inner iterations
                                i_n_branching=par.i_n_branching,  # <-- number of branching subregions
                                i_c=par.i_c,  # <--  increasing number of sampling points ofr all
                                i_replication=par.i_replication,  # <-- initial number of replication
                                i_stopping_max_k=par.i_stopping_max_k,  # <-- maximum number of outer iterations
                                i_max_num_simulation_per_run=par.i_max_num_simulation_per_run,
                                # <-- maximum number of simulations per iteration
                                f_elite_worst_sampling_para=par.f_elite_worst_sampling_para)  # <-- parameters that determine the number of simulation runs for elite and worst subregions


calib_manager = CalibManager(name=name,
                             config_builder=cb,
                             map_sample_to_model_input_fn=map_sample_to_model_input,
                             sites=sites,
                             next_point=optimtool_PBnB, #optimtool,
                             sim_runs_per_param_set=2,
                             max_iterations=20,
                             plotters=plotters)

run_calib_args = {
    "calib_manager": calib_manager
}


if __name__ == "__main__":

    if run_mode == 1:
        # 2. Run these files after input files are set up.  This will not work correctly if the file_setup() line is not commented out, sorry!

        if run_location == "LOCAL":
            SetupParser.init("LOCAL")

        else:
            SetupParser.init()

            SetupParser.set("HPC", "priority", priority)
            SetupParser.set("HPC", "node_group", coreset)

        cm = run_calib_args["calib_manager"]
        cm.run_calibration()

        # comps_exp.submit_experiment(num_seeds=4,simple_intervention_sweep=False)
