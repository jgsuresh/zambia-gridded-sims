"""
Run multi-node simulation
"""

from experiment_setup_new import *
from simtools.SetupParser import SetupParser

# ===================================================================================
exp_name = 'chiyabi_uniform_migrsweep'

catch = "Chiyabi"

gravity_migr_params = np.array([7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00])

base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
grid_pop_csv_file = base + 'data/gridded_pop/cleaned/chiyabi_max_pop.csv'
imm_1node_fp = base + "data/immunity/Sinamalima_1_node_immune_init_p1_33_p2_117.json"

# Intervention file names:
# healthseek_fn = base + 'data/interventions/chiyabi/gridded-uniform/grid_chiyabi_hfca_healthseek_events.csv'
# itn_fn = base + 'data/interventions/chiyabi/gridded-uniform/grid_chiyabi_hfca_itn_events.csv'
# irs_fn = base + 'data/interventions/chiyabi/gridded-uniform/grid_chiyabi_hfca_irs_events.csv'
healthseek_fn = base + 'data/interventions/chiyabi_healthseek.csv'
itn_fn = base + 'data/interventions/chiyabi_itn.csv'
irs_fn = base + 'data/interventions/chiyabi_irs.csv'
msat_fn = base + 'data/interventions/chiyabi/gridded-uniform/grid_chiyabi_hfca_msat_events.csv'
mda_fn = base + 'data/interventions/chiyabi/gridded-uniform/grid_chiyabi_hfca_mda_events.csv'
stepd_fn = base + 'data/interventions/chiyabi/gridded-uniform/grid_chiyabi_hfca_stepd_events.csv'


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
                             rcd_people_num=5,
                             gravity_migr_params=gravity_migr_params,
                             num_cores=6,
                             healthseek_fn=healthseek_fn,
                             itn_fn=itn_fn,
                             irs_fn=irs_fn,
                             msat_fn=msat_fn,
                             mda_fn=mda_fn,
                             stepd_fn=stepd_fn)

# comps_exp.file_setup(generate_climate_files=False) # 1. Run this first to set up input files.  When submitting experiment to COMPS, comment this line out.

if __name__ == "__main__":
    # pass

    # 2. Run these files after input files are set up.  This will not work correctly if the file_setup() line is not commented out, sorry!
    SetupParser.init()
    SetupParser.set("HPC", "priority", "Highest")
    # SetupParser.set("HPC", "node_group", "emod_32cores")
    comps_exp.submit_experiment(num_seeds=5,migration_sweep=True,simple_intervention_sweep=False)
