"""
Run multi-node simulation
"""

from experiment_setup import *


# ===================================================================================
exp_name = 'gravity_test_v0'


base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/data/experiments/{}/inputs/'.format(exp_name)
grid_pop_csv_file = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/data/gridded_pop/cleaned/chiyabi_max_pop.csv'
imm_1node_fp = "C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia-gridded-sims/data/immunity/Sinamalima_1_node_immune_init_p1_33_p2_117.json"

start_year = 2007
num_years = 12

comps_exp = COMPS_Experiment(base,
                             exp_name,
                             grid_pop_csv_file=grid_pop_csv_file,
                             imm_1node_fp=imm_1node_fp,
                             immunity_on=True,
                             migration_on=True,
                             intervention_on=False,
                             gridded_interventions=False,
                             start_year=start_year,
                             sim_length_years=num_years,
                             rcd_people_num=5,
                             # pop_start=p,
                             change_const=False,
                             change_water=True,
                             gravity_migr=True)

# comps_exp.file_setup("SingleNode")#,generate_migration_files=False,generate_demographics_file=True,generate_climate_files=False,generate_immunity_file=True)
# comps_exp.file_setup("MultiNode",generate_migration_files=True)


if __name__ == "__main__":
    # pass
    SetupParser.init()
    SetupParser.set("HPC", "priority", "Normal")
    comps_exp.submit_experiment("MultiNode", num_seeds=2)


    # comps_exp.submit_experiment("SingleNode",num_seeds=10,custom_name='L1_with_summary')
    # comps_exp.submit_experiment("MultiNode",vector_migration_sweep=True)
