"""
Run multi-node simulation of Chiyabi in Lake Kariba region.

"""

from shutil import copyfile
from experiment_setup import *


# ===================================================================================
p=4761
exp_name = 'L10_p{}'.format(p)
base = 'C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/experiments/{}/inputs/'.format(exp_name)
# grid_pop_csv_file = 'C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/gridding/chiyabi-all-rds-max.csv'
grid_pop_csv_file = 'C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/gridding/chiyabi-redo_cait_pop{}.csv'.format(p)
# grid_pop_csv_file = 'C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/gridding/chiyabi-all-rds-max-pop1000.csv'
imm_1node_fp = "C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/inputs/Immunity/Sinamalima_1_node_immune_init_p1_33_p2_117.json"

start_year = 2007
num_years = 12
# start_year = 2001
# num_years = 14


comps_exp = COMPS_Experiment(base, exp_name, grid_pop_csv_file=grid_pop_csv_file, imm_1node_fp=imm_1node_fp,
                             immunity_on=True, migration_on=True, intervention_on=True,
                             start_year=start_year, sim_length_years=num_years,
                             rcd_people_num=5, pop_start=p, change_const=False, change_water=True)

# comps_exp.file_setup("SingleNode")#,generate_migration_files=False,generate_demographics_file=True,generate_climate_files=False,generate_immunity_file=True)
# comps_exp.file_setup("MultiNode")


if __name__ == "__main__":
    pass
    SetupParser.init()
    SetupParser.set("HPC", "priority", "Normal")
    comps_exp.submit_experiment("SingleNode", num_seeds=10,custom_name='L10_p4761_infectivityupdate')
    # comps_exp.submit_experiment("SingleNode",num_seeds=10,custom_name='L1_with_summary')

    # comps_exp.submit_experiment("MultiNode",vector_migration_sweep=True)


# L10: constant fixed, water_vegetation scales with population
