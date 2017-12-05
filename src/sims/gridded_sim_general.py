import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt

# Place to keep just generally useful functions



def find_cells_for_this_catchment(catch, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    # Find which grid cells correspond to a given HFCA
    df = pd.read_csv(base + 'data/interventions/kariba/2017-11-27/raw/grid_lookup.csv')

    if catch == 'all':
        return np.array(df['grid_cell'])
    else:
        df_catch = df[df['catchment']==catch.lower()]
        return np.array(df_catch['grid_cell'])

def convert_from_dtk_node_ids_to_grid_cells_using_demo(dtk_node_ids,demo):
    # Parse a demographics file to return grid cell IDs that correspond to given simulation node IDs:
    with open(demo,'r') as f:
        demo_df = json.load(f)

    nnodes_demo = len(demo_df['Nodes'])

    node_ids = []
    grid_cell_ids = []
    for node in demo_df['Nodes']:
        node_ids.append(node['NodeID'])
        grid_cell_ids.append(node['FacilityName'])

    # node_ids = pd.Series(node_ids)
    # grid_cell_ids = pd.Series(grid_cell_ids)

    demo_df = pd.DataFrame({
        'node_ids': node_ids,
        'grid_cell_ids': grid_cell_ids
    })

    given_node_df = pd.DataFrame({
        'dtk_node_ids': dtk_node_ids
    })

    # Use merge to get the grid cell IDs that correspond to the requested dtk node IDs, then return as an array
    full_df = given_node_df.merge(demo_df,how='left',left_on='dtk_node_ids',right_on='node_ids')
    return np.array(full_df['grid_cell_ids'])


def get_lat_long_grid_cells(cell_ids, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    df = pd.read_csv(base + 'data/interventions/kariba/2017-11-27/raw/grid_lookup.csv')

    lat = search_dataframe(df,"grid_cell",cell_ids,"mid_y")
    lon = search_dataframe(df,"grid_cell",cell_ids,"mid_x")

    return [lat,lon]

def get_lat_long_dtk_nodes(dtk_node_ids, demo, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    cell_ids = convert_from_dtk_node_ids_to_grid_cells_using_demo(dtk_node_ids,demo)
    [lat,lon] = get_lat_long_grid_cells(cell_ids, base=base)
    return [lat,lon]

def search_dataframe(df,search_col,search_vals,return_col):
    # Return values from "return_col" which are from corresponding rows of "search_vals" in "search_col", IN ORDER

    # Merge the two:
    temp_df = pd.DataFrame({
        search_col: search_vals
    })

    search_df = temp_df.merge(df,how='left',left_on=search_col,right_on=search_col)
    return np.array(search_df[return_col])


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
    fn = base + "data/immunity/{}_1_node/{}_1_node_immune_init_p1_{}_p2_{}.json".format(climate_category,climate_category,nearest_funest_param,nearest_arab_param)
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


def scatter_lat_long_on_map(lon,lat,
                            C=None,S=None,cbar_label=None, cmap=False, clim=None,
                            lat_range=None, lon_range=None,savefig=False,title=None):
    # lat, lon are coordinates of data
    # C is some third data quantity used to color points
    # S is some fourth quantity used to size points
    import mpl_toolkits.basemap as base

    # fig, ax = plt.subplots()
    fig = plt.figure(figsize=(10,10))
    ax=plt.subplot()

    if isinstance(lat_range,list) and isinstance(lon_range,list):
        x_min = lon_range[0]
        x_max = lon_range[1]
        y_min = lat_range[0]
        y_max = lat_range[1]
    else:
        x_min = np.min(lon)
        x_max = np.max(lon)
        y_min = np.min(lat)
        y_max = np.max(lat)

        w, h = x_max - x_min, y_max - y_min

        x_min = x_min - 0.1*w
        x_max = x_max + 0.1*w
        y_min = y_min - 0.1*h
        y_max = y_max + 0.1*h

    # Get/plot the background map:
    m = base.Basemap(
        # projection = 'merc',
        # ellps = 'WGS84',
        llcrnrlon=x_min,
        llcrnrlat=y_min,
        urcrnrlon=x_max,
        urcrnrlat=y_max,
        lat_ts=0,
        epsg=4269
    )

    # see http://server.arcgisonline.com/arcgis/rest/services (section services) for more options
    m.arcgisimage(service='ESRI_Imagery_World_2D', xpixels=700, verbose=True)
    # m.arcgisimage(service='NatGeo_World_Map', xpixels = 2000, verbose= True)
    # m.drawcountries(linewidth=.25,linestyle='solid')

    min_x, min_y = m(x_min, y_min)
    max_x, max_y = m(x_max, y_max)
    corr_w, corr_h = max_x - min_x, max_y - min_y

    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)

    ax.set_aspect(1)

    # Scatter plot data:
    if isinstance(S,int):
        S = np.ones_like(lon) * S
    elif not isinstance(S,np.ndarray):
        S = np.ones_like(lon)*5

    if not isinstance(C,np.ndarray):
        C = 'C0'

    if not cmap:
        cmap = plt.cm.viridis

    if not isinstance(clim,list):
        clim = [min(C),max(C)]

    sc = ax.scatter(lon, lat, c=C, vmin=clim[0], vmax=clim[1],marker='s', s=S, cmap=cmap,edgecolors='black')

    plt.colorbar(sc)
    plt.title(title)

    if savefig == False:
        plt.show()
    else:
        plt.savefig(savefig)

    return ax