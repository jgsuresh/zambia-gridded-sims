import pandas as pd
import numpy as np


def find_milen_cluster_for_grid_cells(cells):
    # Load Caitlin's cluster_to_grid_lookup.csv to identify, for each grid cell, which Milen-cluster it belongs to
    base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
    cl_lookup = pd.read_csv(base + "data/interventions/kariba/2017-11-27/raw/cluster_to_grid_lookup.csv")

    #


    pass
# Load Milen's best-fit larval params JSON file
# Given

# larval_params_for_grid_cells