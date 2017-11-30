import pandas as pd
import numpy as np
import json

# Place to keep just generally useful functions



def find_cells_for_this_catchment(catch, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    # Find which grid cells correspond to a given HFCA
    df = pd.read_csv(base + 'data/interventions/kariba/2017-11-27/raw/grid_lookup.csv')

    if catch == 'all':
        return np.array(df['grid_cell'])
    else:
        df_catch = df[df['catchment']==catch.lower()]
        return np.array(df_catch['grid_cell'])


#########################################################################################
# Functions relating to connections between grid cells and Milen's clusters:

def find_milen_cluster_for_grid_cells(cell_ids, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    # Load Caitlin's cluster_to_grid_lookup.csv to identify, for each grid cell, which Milen-cluster it belongs to
    #NOTE: Assumes that cells is sorted
    cl_lookup = pd.read_csv(base + "data/milen_clusters/cluster_to_grid_lookup.csv")

    # fixme This assumes that every possible cell in cell_ids is in the lookup_cells.  Fails if not the case.
    # lookup_cells = np.array(cl_lookup['grid_cell'])
    # lookup_id = np.array(cl_lookup['cluster_id'])
    # in_desired_cells = np.in1d(lookup_cells,cell_ids)
    # cluster_ids = list(lookup_id[in_desired_cells])

    hold_df = pd.DataFrame({
        "cell_ids": cell_ids
    })

    hold_df = hold_df.merge(cl_lookup,how='left',left_on='cell_ids',right_on='grid_cell')
    cluster_ids = hold_df['cluster_id']
    cluster_ids = cluster_ids.fillna(method='ffill') # Interpolate for any missing cells
    cluster_ids = cluster_ids.fillna(method='bfill')
    cluster_ids = list(cluster_ids)

    # print cell_ids
    # print cluster_ids
    return cluster_ids

def milen_cluster_larval_param_fit(cluster_id, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    # Load Milen's best-fit larval params JSON file and get larval param fit for this cluster_id
    fn = base + "data/larval_params/milen_best_fits.json"
    f = open(fn,"r")
    larval_fits_dict = json.load(f)
    f.close()

    return_dict = larval_fits_dict[cluster_id]['fit']['params'].copy()
    del return_dict['drug_cov']
    del return_dict['itn_level']

    return return_dict

def find_milen_larval_param_fit_for_grid_cells(cell_ids):
    # Return the best-fit larval parameters for a set of grid cell IDs
    cluster_ids = find_milen_cluster_for_grid_cells(cell_ids)

    arab_params = np.zeros_like(cluster_ids)
    funest_params = np.zeros_like(cluster_ids)

    i = 0
    for cl_id in cluster_ids:
        # print "cl_id ",cl_id
        param_dict = milen_cluster_larval_param_fit(cl_id)
        arab_params[i] = param_dict['arabiensis_sc']
        funest_params[i] = param_dict['funestus_sc']

        i += 1

    arab_params = arab_params.astype(np.float)
    funest_params = funest_params.astype(np.float)
    return [arab_params,funest_params]

def milen_cluster_climate_category(cluster_id, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    # Return the climate category for a given cluster ID
    lookup = pd.read_csv(base + "data/milen_clusters/cluster_climate_mapping.csv")
    return lookup[lookup['cluster']==cluster_id]['climate.cat'].item()

def closest_milen_immunity_overlay_filename(arab_param,funest_param,climate_category, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    # Find filename of immunity file which is closest match, within given climate category, to specified arabiensis/funestus larval habitat params
    def find_nearest_from_precomputed(x,precomputed_x):
        precomp_x_arr = np.array(precomputed_x) #assume precomputed_x is a list
        closest = np.argmin(np.abs(x-precomp_x_arr))
        return precomputed_x[closest]

    # Get arab and funest precomputed values for specified climate category
    precomputed_arab_params, precomputed_funest_params = get_arab_funest_precomputed_values_for_climate_category(climate_category)

    nearest_arab_param = find_nearest_from_precomputed(arab_param,precomputed_arab_params)
    nearest_funest_param = find_nearest_from_precomputed(funest_param, precomputed_funest_params)

    # Return filename for these parameters:
    fn = base + "data/immunity/{}_1_node/{}_1_node_immune_init_p1_{}_p2_{}.json".format(climate_category,climate_category,nearest_arab_param,nearest_funest_param)
    return fn

def get_arab_funest_precomputed_values_for_climate_category(climate_category, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    # Get all filenames in that climate folder, and extract the precomputed vector param values by parsing the filenames

    def represents_int(s):
        try:
            int(s)
            return True
        except ValueError:
            return False

    def get_arab_funest_from_filename(fn):
        # Extract the two parameters from the file name:
        break1 = fn.split('_p1_')
        break1 = break1[1]

        break2 = break1.split('_p2_')
        if represents_int(break2[0]):
            funest_param = int(break2[0])
        else:
            funest_param = np.float(break2[0])

        break2 = break2[1]
        break3 = break2.split('.json')
        if represents_int(break3[0]):
            arab_param = int(break3[0])
        else:
            arab_param = np.float(break3[0])

        return [arab_param, funest_param]

    import glob
    folder_name = base + "data/immunity/{}_1_node".format(climate_category)
    file_list = glob.glob(folder_name + "/*.json")

    arab_params = [] #np.zeros(len(file_list))
    funest_params = [] #np.zeros(len(file_list))

    # Loop over all filenames and extract vector parameters from the filenames
    # i = 0
    # print file_list
    for fn in file_list:
        a,f = get_arab_funest_from_filename(fn)
        # arab_params[i] = a
        # funest_params[i] = f
        arab_params.append(a)
        funest_params.append(f)

    # return [np.unique(arab_params),np.unique(funest_params)]
    unique_arab_params = list(set(arab_params))
    unique_funest_params = list(set(funest_params))

    unique_arab_params.sort()
    unique_funest_params.sort()
    return [unique_arab_params,unique_funest_params]

def closest_milen_immunity_overlay_filenames_for_grid_cells(cell_ids):
    cluster_ids = find_milen_cluster_for_grid_cells(cell_ids)

    fn_list = []
    for cluster_id in cluster_ids:
        larval_param_dict = milen_cluster_larval_param_fit(cluster_id)
        climate_category = milen_cluster_climate_category(cluster_id)

        fn = closest_milen_immunity_overlay_filename(larval_param_dict['arabiensis_sc'],
                                                     larval_param_dict['funestus_sc'],
                                                     climate_category)

        fn_list.append(fn)

    return fn_list




############################################################################################################
# def generate_immun_dict_from_demo_file(cell_ids, node_ids):
#     n_nodes = len(cell_ids)
#     immun_fn_list = closest_milen_immunity_overlay_filenames_for_grid_cells(cell_ids)
#
#     d = {}
#     d["Nodes"] = []
#     d["Metadata"] = {}
#     d["Metadata"]["Author"] = "Josh Suresh"
#     d["Metadata"]["IdReference"] = "Gridded world grump30arcsec"
#     d["Metadata"]["NodeCount"] = n_nodes
#
#     for i in range(n_nodes):
#         immun_fn = immun_fn_list[i]
#         f = open(immun_fn,'r')
#         imm_dict = json.load(f)
#         f.close()
#
#         node = {}
#         node["NodeAttributes"] = imm_dict['Defaults'].copy()
#         node["NodeID"] = node_ids[i]
#         d["Nodes"].append(node)

