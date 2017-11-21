import pandas as pd
import numpy as np

# Open one of Caitlin's intervention files.
# Get grid cell, and catchment name.  save this as an intermediate file.
# For a given catchment name, get grid cells which correspond to this,
# And generate IRS, ITN, and healthseek files for that particular catchment.


def get_grid_cells_for_catch(catch):
    # Open Caitlin's intervention files and return which grid cells correspond to a given catchment name
    base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
    fn1 = base + 'data/interventions/grid_all_irs_events.csv'
    fn2 = base + 'data/interventions/grid_all_itn_events.csv'

    def return_grid_cells_this_catch(fn):
        df = pd.read_csv(fn)
        if catch == 'all':
            raise Exception(NotImplementedError)
        else:
            df_catch = df[df['catch'] == catch]

        grid_cell_numbers = np.unique(df_catch['grid_cell'])
        return grid_cell_numbers

    grid_nums_1 = return_grid_cells_this_catch(fn1)
    grid_nums_2 = return_grid_cells_this_catch(fn2)

    grid_cell_numbers = np.union1d(grid_nums_1,grid_nums_2)
    return grid_cell_numbers

def get_grid_cells_for_chiyabi():
    # Open Caitlin's intervention files and return which grid cells correspond to a given catchment name
    base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
    fn = base + 'data/interventions/chiyabi/gridded-uniform/grid_chiyabi_hfca_irs_events.csv'

    df = pd.read_csv(fn)
    grid_cell_numbers = np.unique(df['grid.cell'])
    return grid_cell_numbers


def generate_intervention_files_by_catchment(catch):
    base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'

    # grid_cell_nums = get_grid_cells_for_catch(catch)
    grid_cell_nums = get_grid_cells_for_chiyabi()

    def create_subfile_for_catch(fn,type):
        df = pd.read_csv(fn)
        df_catch = df[np.in1d(df['grid_cell'],grid_cell_nums)]
        df_catch = df_catch.reset_index()
        df_catch = df_catch.drop('index',axis=1)
        df_catch.rename(columns={'grid_cell':'grid.cell'},inplace=True)

        # Save file
        outfile = base + 'data/interventions/{}_{}.csv'.format(catch.lower(),type)
        df_catch.to_csv(outfile)

    create_subfile_for_catch(base + 'data/interventions/grid_all_healthseek_events.csv','healthseek')
    create_subfile_for_catch(base + 'data/interventions/grid_all_irs_events.csv', 'irs')
    create_subfile_for_catch(base + 'data/interventions/grid_all_itn_events.csv', 'itn')

if __name__ == "__main__":
    generate_intervention_files_by_catchment('Chiyabi')