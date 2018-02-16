
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

from dtk.vector.study_sites import configure_site

from GriddedCalibSite import GriddedCalibSite

import matplotlib
matplotlib.use('Agg')


# Workflow to-do (Sun morning, 2/11)
# Generate some serialized files from the calib-week calib script.  Put these in start_files/
# Make some fake dict hfca_clusters.json, put it in the folder
# Make sure param-generating code is correct
# Make sure map_sample_to_model_input is correct (larval params only on particular nodeset)
#

# Calibration parameters
calib_folder_name='TestSerialize'

# This folder must contain serialized files in start_files/ folder, and hfca_clusters.json dictionary file which groups nodes into clusters

run_mode = 0
serialize = True

#============================================================================================================
# Run parameters:
catch = "Munyumbwe"
exp_name = '{}_calibtest'.format(catch)

base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
grid_pop_csv_file = base + 'data/mozambique/grid_population.csv' # Max pop in grid cell for entire region

# Intervention file names:
healthseek_fn = base + 'data/mozambique/grid_all_healthseek_events.csv'
itn_fn = base + 'data/mozambique/grid_all_itn_events.csv'
irs_fn = base + 'data/mozambique/grid_all_irs_events.csv'
# msat_fn = base + 'data/mozambique/grid_all_msat_events.csv'
mda_fn = base + 'data/mozambique/grid_all_mda_events.csv'
# stepd_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_stepd_events.csv'

start_year = 1994
num_years = 24
serialize_year = 2014

num_cores = 8

#============================================================================================================
# Set up run:

gravity_migr_params = np.array([7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00])

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
    comps_exp.file_setup()  # 1. Run this first to set up input files.

cb = comps_exp.return_cb_for_calibration()

if serialize:
    cb.update_params({'Serialization_Time_Steps': [365 * (serialize_year-start_year)]})

#============================================================================================================
# Set up calibration:
# Calibration-specific stuff:
sites = [GriddedCalibSite()]

plotters = [LikelihoodPlotter(combine_sites=True),
            SiteDataPlotter(num_to_plot=5, combine_sites=True), # necessary because it generates the LL_all.csv file
            OptimToolPlotter()  # OTP must be last because it calls gc.collect()
            ]

#fixme Load hfca_clusters.json to find # of clusters and which nodes are in which cluster:

# for cluster in list_of_clusters:
#     add params to this list:
params = [
    {
        'Name': 'funestus_scale',
        'Dynamic': True,
        'MapTo': 'funestus_scale', # <-- DEMO: Custom mapping, see map_sample_to_model_input below
        'Guess': 5,
        'Min': 0.1,
        'Max': 50
    },
    {
        'Name': 'arabiensis_scale',
        'Dynamic': True,
        'MapTo': 'arabiensis_scale',
        'Guess': 35,
        'Min': 0.1,
        'Max': 100
    },
]


def map_sample_to_model_input(cb, sample):
    tags = {}
    sample = copy.deepcopy(sample)

    if 'arabiensis_scale' in sample:
        a_sc = sample.pop('arabiensis_scale')

        hab = {'arabiensis': {'TEMPORARY_RAINFALL': 1e8 * a_sc, 'CONSTANT': 2e6}}
        set_larval_habitat(cb, hab)

        tags.update({'arabiensis_scale': a_sc})

    if 'funestus_scale' in sample:
        f_sc = sample.pop('funestus_scale')

        hab = {'funestus': {
            "WATER_VEGETATION": 2e6, # Caitlin.  Milen had 2e7,
            "LINEAR_SPLINE": {
                "Capacity_Distribution_Per_Year": {
                    "Times": [0.0, 30.417, 60.833, 91.25, 121.667, 152.083, 182.5, 212.917, 243.333, 273.75, 304.167,
                              334.583],
                    "Values": [0.0, 0.0, 0.0, 0.2, 0.8, 1.0, 1.0, 1.0, 0.5, 0.2, 0.0, 0.0] # Caitlin
                    # "Values": [0.0, 0.0, 0.0, 0.0, 0.2, 1.0, 1.0, 1.0, 0.5, 0.2, 0.0, 0.0] # Milen
                },
                "Max_Larval_Capacity": 1e8 * f_sc
            }
        }
        }

        set_larval_habitat(cb, hab)

        tags.update({'funestus_scale': f_sc})

    for name,value in sample.items():
        print('UNUSED PARAMETER:', name)
    assert( len(sample) == 0 ) # All params used

    #fixme Draw from serialized files at specific time:
    # If you can find LL_all, then use that to find best fit so far, and use its serialized files.
    # If not, then use the serialized files in start_files instead.


    return tags


optimtool = OptimTool(params,
                      samples_per_iteration=5,
                      center_repeats=1)

calib_manager = CalibManager(name='TestOptimMotaze',
                             config_builder=cb,
                             map_sample_to_model_input_fn=map_sample_to_model_input,
                             sites=sites,
                             next_point=optimtool,
                             sim_runs_per_param_set=2,
                             max_iterations=8,
                             plotters=plotters)

run_calib_args = {
    "calib_manager": calib_manager
}
