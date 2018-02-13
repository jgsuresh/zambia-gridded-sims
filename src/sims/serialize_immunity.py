
import os
import json
import pandas as pd

from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.utils.reports.VectorReport import add_human_migration_report

from simtools.Utilities.COMPSUtilities import COMPS_login
from simtools.SetupParser import SetupParser
from simtools.Utilities.COMPSUtilities import translate_COMPS_path
from simtools.DataAccess.ExperimentDataStore import ExperimentDataStore
from simtools.ModBuilder import ModBuilder, ModFn

from dtk.vector.species import set_params_by_species, update_species_param

from dtk.interventions.health_seeking import add_health_seeking
from dtk.interventions.itn_age_season import add_ITN_age_season

from malaria.interventions.malaria_drug_campaigns import add_drug_campaign
from malaria.reports.MalariaReport import add_event_counter_report


# Write up #1
# Assume: starting from a completely new place that we have no "reasonable guess for".
# If we have no reasonable guess:
    # Do a simple 2D calibration everywhere, with 20 year immunity burn-ins.  Three iterations.  (e.g. for Munyumbwe, this would be about 6 hours)
    # Serialize all files pre-intervention.
    # Use the best run from this as the starting point for what follows
# Elif we do have a reasonable guess:
    # Run a long burn-in simulation with those parameters.
    # Serialize this run pre-intervention.  Use this as the starting point for what follows

# Use the serialized pre-intervention files to start a new calibration with all parameters of interest (e.g. multiple parameters for different parts of the settlement)
# that has a very short burn-in (~5 years).
# Set up the calibration so that ALL of these runs will serialize pre-intervention



# Write up #2:
# Load serialized files from "good enough guess"
serialization_date = 20*365
serialization_exp_id = "95c9bd38-72ae-e711-9414-f0921c16b9e5"
num_cores = 8

serialize = True  # If true, save serialized files for calibration runs


COMPS_login("https://comps.idmod.org")
expt = ExperimentDataStore.get_most_recent_experiment(serialization_exp_id)

df = pd.DataFrame([x.tags for x in expt.simulations])
df['outpath'] = pd.Series([sim.get_path() for sim in expt.simulations])

builder = ModBuilder.from_list([[
    ModFn(DTKConfigBuilder.set_param, 'Serialized_Population_Path', '{path}/output'.format(path=df['outpath'][x])),
    ModFn(DTKConfigBuilder.set_param, 'Run_Number', df['Run_Number'][x])
                                ]
    for x in df.index

])
if num_cores > 1:
    cb.update_params({'Serialized_Population_Filenames':
                          ['state-%05d-%03d.dtk' % (serialization_date, corenum) for corenum in
                           range(num_cores)]
                      })
else:
    cb.update_params({'Serialized_Population_Filenames':
                          ['state-%05d.dtk' % serialization_date]
                      })



# Load in file that tells us what nodes correspond to what subclusters
# f is a dictionary whose keys are node IDs, and whose values are subcluster IDs.  Header has N_subclusters
# f = pd.read_json("hfca_clusters.json")

# Create custom mapping function based on how many subclusters we have.  Will be 2*N_subcl parameters to optimize over.

# In each optimization, start from serialized files + 5 years before interventions begin.
# Each optimization run also serializes its own population pre-intervention.

# Continue this until Calibration is done.