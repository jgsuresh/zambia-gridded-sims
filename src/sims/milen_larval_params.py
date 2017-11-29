import pandas as pd
import numpy as np
import json


def find_milen_cluster_for_grid_cells(cells):
    # Load Caitlin's cluster_to_grid_lookup.csv to identify, for each grid cell, which Milen-cluster it belongs to
    #fixme Assumes that cells is sorted
    base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
    cl_lookup = pd.read_csv(base + "data/interventions/kariba/2017-11-27/raw/cluster_to_grid_lookup.csv")

    lookup_cells = np.array(cl_lookup['grid_cell'])
    lookup_id = np.array(cl_lookup['cluster_id'])
    in_desired_cells = np.in1d(lookup_cells,cells)

    cluster_ids = list(lookup_id[in_desired_cells])

    return cluster_ids

def milen_larval_param_fit_for_cluster(cluster_id):
    # Load Milen's best-fit larval params JSON file and get larval param fit for this cluster_id
    base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
    fn = base + "data/larval_params/milen_best_fits.json"
    f = open(fn,"r")
    larval_fits_dict = json.load(f)
    f.close()

    return_dict = larval_fits_dict[cluster_id]['fit']['params']
    del return_dict['drug_cov']
    del return_dict['itn_level']

    return return_dict

def find_milen_larval_param_fit_for_grid_cells(cells):
    cluster_ids = find_milen_cluster_for_grid_cells(cells)

    arab_params = np.zeros_like(cluster_ids)
    funest_params = np.zeros_like(cluster_ids)

    i = 0
    for cl_id in cluster_ids:
        param_dict = milen_larval_param_fit_for_cluster(cl_id)
        arab_params[i] = param_dict['arabiensis_sc']
        funest_params[i] = param_dict['funestus_sc']

    return [arab_params,funest_params]



# if __name__ == "__main__":
#     base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
#     cl_lookup = pd.read_csv(base + "data/interventions/kariba/2017-11-27/raw/cluster_to_grid_lookup.csv")
