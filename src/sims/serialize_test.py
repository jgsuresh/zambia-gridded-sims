
"""
Test serialization workflow.  Starts a calibrator which generates serialized files, then runs itself from the best version of these serialized files
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

from dtk.vector.study_sites import configure_site

from GriddedCalibSite import GriddedCalibSite

import matplotlib
matplotlib.use('Agg')

run_mode = 1
run_location = "HPC" #""LOCAL"
priority = "AboveNormal"
coreset = "emod_abcd" #"emod_32cores"
num_cores = 8
milen_catch_list = ['bbondo', 'chabbobboma', 'chisanga', 'chiyabi', 'luumbo', 'munyumbwe', 'nyanga chaamwe',
                    'sinafala', 'sinamalima']
catch = milen_catch_list[5]


# ===================================================================================

exp_name = '{}_calibtest_2'.format(catch)

gravity_migr_params = np.array([7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00])

base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
grid_pop_csv_file = base + 'data/gridded_pop/cleaned/all_max_pop.csv' # Max pop in grid cell for entire region

# Intervention file names:
healthseek_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_healthseek_events.csv'
itn_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_itn_events.csv'
irs_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_irs_events.csv'
msat_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_msat_events.csv'
mda_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_mda_events.csv'
stepd_fn = base + 'data/interventions/kariba/2017-11-27/raw/grid_all_stepd_events.csv'


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
                             rcd_people_num=10,
                             gravity_migr_params=gravity_migr_params,
                             num_cores=num_cores,
                             healthseek_fn=healthseek_fn,
                             itn_fn=itn_fn,
                             irs_fn=irs_fn,
                             msat_fn=msat_fn,
                             mda_fn=mda_fn,
                             stepd_fn=stepd_fn,
                             larval_params_mode="calibrate",
                             immunity_mode="naive")

if run_mode == 0:
    comps_exp.file_setup()  # 1. Run this first to set up input files.  When submitting experiment to COMPS, comment this line out.

cb = comps_exp.return_cb_for_calibration()


#fixme EVENTUALLY, ADD CLUSTER-CHUNKING ALGORITHM HERE.  FOR NOW, CALIBRATE OVER ALL CELLS
# Break the file into clusters
# Loop over clusters.  For each cluster, calibrate:




# Calibration-specific stuff:
sites = [GriddedCalibSite()]
# The default plotters used in an Optimization with OptimTool
plotters = [LikelihoodPlotter(combine_sites=True),
            SiteDataPlotter(num_to_plot=5, combine_sites=True), #needed for LL_all.csv
            OptimToolPlotter()  # OTP must be last because it calls gc.collect()
            ]

params = [
    {
        'Name': 'funestus_scale',
        'Dynamic': True,
        'MapTo': 'funestus_scale', # <-- DEMO: Custom mapping, see map_sample_to_model_input below
        'Guess': 25,
        'Min': 0.1,
        'Max': 150
    },
    {
        'Name': 'arabiensis_scale',
        'Dynamic': True,
        'MapTo': 'arabiensis_scale',
        'Guess': 25,
        'Min': 0.1,
        'Max': 150
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


    # Search for LL_all.csv.  If it does not exist, then this is the first run.  If it does, then this is a later run.
    if os.path.isfile(base + 'src/sims/TestSerialize/_plots/LL_all.csv'):
        cb.update_params({'Serialization_Time_Steps': [365 * 5]})
    else:
        cb.update_params({'Serialization_Time_Steps': [365 * 20]})
        # cb.update_params({'Serialized_Population_Filenames':
        #                       ['state-%05d-%03d.dtk' % (serialization_date, corenum) for corenum in
        #                        range(num_cores)]
        #                   })



    for name,value in sample.items():
        print('UNUSED PARAMETER:', name)
    assert( len(sample) == 0 ) # All params used

    return tags



optimtool = OptimTool(params,
                      samples_per_iteration=25,
                      center_repeats=1)

calib_manager = CalibManager(name='TestOptimSerialize',
                             config_builder=cb,
                             map_sample_to_model_input_fn=map_sample_to_model_input,
                             sites=sites,
                             next_point=optimtool,
                             sim_runs_per_param_set=1,
                             max_iterations=8,
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
