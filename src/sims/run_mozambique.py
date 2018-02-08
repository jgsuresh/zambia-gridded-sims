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

from dtk.vector.study_sites import configure_site

from GriddedCalibSite import GriddedCalibSite

import matplotlib
matplotlib.use('Agg')


run_mode = 1
run_location = "HPC" #""LOCAL"
priority = "AboveNormal"
coreset = "emod_abcd" #"emod_32cores"
num_cores = 8
catch = "Motaze"

# ===================================================================================

exp_name = '{}_seeds'.format(catch)

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

larval_params = {
    "const_h": 1,
    "temp_h": 35,
    "water_h": 1,
    "linear_h": 1
}


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
                             larval_params_mode="uniform",
                             immunity_mode="naive")

if run_mode == 0:
    comps_exp.file_setup(larval_params=larval_params)  # 1. Run this first to set up input files.  When submitting experiment to COMPS, comment this line out.

cb = comps_exp.return_cb_for_calibration()

#fixme EVENTUALLY, ADD CLUSTER-CHUNKING ALGORITHM HERE.  FOR NOW, CALIBRATE OVER ALL CELLS
# Break the file into clusters
# Loop over clusters.  For each cluster, calibrate:



if __name__ == "__main__":

    if run_mode == 1:
        # 2. Run these files after input files are set up.  This will not work correctly if the file_setup() line is not commented out, sorry!

        if run_location == "LOCAL":
            SetupParser.init("LOCAL")

        else:
            SetupParser.init()

            SetupParser.set("HPC", "priority", priority)
            SetupParser.set("HPC", "node_group", coreset)

        comps_exp.submit_experiment(num_seeds=8, simple_intervention_sweep=False)

