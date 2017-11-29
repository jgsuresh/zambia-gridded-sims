import pandas as pd
import numpy as np

# Place to keep just generally useful functions



def find_catchment_cells(catch, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    df = pd.read_csv(base + 'data/interventions/kariba/2017-11-27/raw/grid_lookup.csv')

    if catch == 'all':
        return np.array(df['grid_cell'])
    else:
        df_catch = df[df['catchment']==catch.lower()]
        return np.array(df_catch['grid_cell'])

