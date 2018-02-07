import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt

from relative_time import *

# Place to keep just generally useful functions



def find_cells_for_this_catchment(catch,
                                  base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/',
                                  path_from_base='data/interventions/kariba/2017-11-27/raw/grid_lookup.csv'):
    # Find which grid cells correspond to a given HFCA
    df = pd.read_csv(base + path_from_base)

    if catch == 'all':
        return np.array(df['grid_cell'])
    else:
        df_catch = df[df['catchment'] == catch]
        # df_catch = df[np.logical_or(df['catchment'] == catch.capitalize(),
        #                             df['catchment'] == catch.lower())]
        return np.array(df_catch['grid_cell'])


def parse_demo_file(demo):
    if type(demo) is str:
        with open(demo, 'r') as f:
            demo_df = json.load(f)
    elif type(demo) is dict:
        demo_df = demo

    node_ids = []
    grid_cell_ids = []
    for node in demo_df['Nodes']:
        node_ids.append(node['NodeID'])
        grid_cell_ids.append(node['NodeAttributes']['FacilityName'])

    # Convert grid_cell_ids to array of integers:
    grid_cell_ids = np.array(grid_cell_ids, dtype=int)

    demo_df = pd.DataFrame({
        'node_ids': node_ids,
        'grid_cell_ids': grid_cell_ids
    })

    return demo_df


def convert_from_dtk_node_ids_to_grid_cells_using_demo(dtk_node_ids, demo):
    # Parse a demographics file to return grid cell IDs that correspond to given simulation node IDs:
    demo_df = parse_demo_file(demo)

    given_node_df = pd.DataFrame({
        'dtk_node_ids': dtk_node_ids
    })

    # Use merge to get the grid cell IDs that correspond to the requested dtk node IDs, then return as an array
    full_df = given_node_df.merge(demo_df, how='left', left_on='dtk_node_ids', right_on='node_ids')
    return np.array(full_df['grid_cell_ids'])


def convert_from_grid_cells_to_dtk_node_ids_using_demo(grid_cells, demo):
    # Parse a demographics file to return grid cell IDs that correspond to given simulation node IDs:
    demo_df = parse_demo_file(demo)

    given_cell_df = pd.DataFrame({
        'grid_cells': grid_cells
    })

    # Use merge to get the grid cell IDs that correspond to the requested dtk node IDs, then return as an array
    full_df = given_cell_df.merge(demo_df, how='left', left_on='grid_cells', right_on='grid_cell_ids')
    return np.array(full_df['node_ids'])


def get_lat_long_grid_cells(cell_ids, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    df = pd.read_csv(base + 'data/interventions/kariba/2017-11-27/raw/grid_lookup.csv')

    lat = search_dataframe(df, "grid_cell", cell_ids, "mid_y")
    lon = search_dataframe(df, "grid_cell", cell_ids, "mid_x")

    return [lat, lon]


def get_lat_long_dtk_nodes(dtk_node_ids, demo, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    cell_ids = convert_from_dtk_node_ids_to_grid_cells_using_demo(dtk_node_ids, demo)
    # print "cell_ids ",cell_ids
    [lat, lon] = get_lat_long_grid_cells(cell_ids, base=base)
    return [lat, lon]


def search_dataframe(df, search_col, search_vals, return_col):
    # Return values from "return_col" which are from corresponding rows of "search_vals" in "search_col", IN ORDER

    # Merge the two:
    temp_df = pd.DataFrame({
        search_col: search_vals
    })

    search_df = temp_df.merge(df, how='left', left_on=search_col, right_on=search_col)
    return np.array(search_df[return_col])


def get_RDT_ref_data_for_grid_cells(grid_cells,
                                    format="combine",
                                    base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/',
                                    path_from_base="data/prevalence/2018-01-23/raw/grid_prevalence_with_dates.csv"):
    # Return grid_cell ID, date, population, and prevalence for all grid cells
    # If format == "combine", return these all as a single list

    # Open relevant file:
    prev_df = pd.read_csv(base + path_from_base)

    in_cells = np.in1d(prev_df["grid_cell"],grid_cells)
    return prev_df[in_cells]

#########################################################################################
# Functions relating to connections between grid cells and Milen's clusters:

def find_milen_cluster_for_grid_cells(cell_ids, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    # Load Caitlin's cluster_to_grid_lookup.csv to identify, for each grid cell, which Milen-cluster it belongs to
    # NOTE: Assumes that cells is sorted
    cl_lookup = pd.read_csv(base + "data/milen_clusters/cluster_to_grid_lookup.csv")

    # fixme This assumes that every possible cell in cell_ids is in the lookup_cells.  Fails if not the case.
    # lookup_cells = np.array(cl_lookup['grid_cell'])
    # lookup_id = np.array(cl_lookup['cluster_id'])
    # in_desired_cells = np.in1d(lookup_cells,cell_ids)
    # cluster_ids = list(lookup_id[in_desired_cells])

    hold_df = pd.DataFrame({
        "cell_ids": cell_ids
    })

    hold_df = hold_df.merge(cl_lookup, how='left', left_on='cell_ids', right_on='grid_cell')
    cluster_ids = hold_df['cluster_id']
    cluster_ids = cluster_ids.fillna(method='ffill')  # Interpolate for any missing cells
    cluster_ids = cluster_ids.fillna(method='bfill')
    cluster_ids = list(cluster_ids)

    # print cell_ids
    # print cluster_ids
    return cluster_ids


def milen_cluster_larval_param_fit(cluster_id, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    # Load Milen's best-fit larval params JSON file and get larval param fit for this cluster_id
    fn = base + "data/larval_params/milen_best_fits.json"
    f = open(fn, "r")
    larval_fits_dict = json.load(f)
    f.close()

    return_dict = larval_fits_dict[cluster_id]['fit']['params'].copy()
    del return_dict['drug_cov']
    del return_dict['itn_level']

    return return_dict


def find_milen_larval_param_fit_for_grid_cells(cell_ids, fudge_milen_habitats=False):
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

        if fudge_milen_habitats:
            [arab_params[i],funest_params[i]]= fudge_milen_habitats_by_hand(cl_id,
                                                                             arab_params[i].astype(np.float),
                                                                             funest_params[i].astype(np.float))

        i += 1

    arab_params = arab_params.astype(np.float)
    funest_params = funest_params.astype(np.float)
    return [arab_params, funest_params]


def fudge_milen_habitats_by_hand(cl_id,arab_param, funest_param):
    # Modify larval habitat parameters by hand
    mod_factor_dict = {
        # Bbondo:
        "80201_1": 0.5,
        "80201_2": 0.5,
        "80201_5": 0.5,
        "80201_6": 0.5,
        "80201_8": 0.5,
        "80201_10": 0.5,
        "80201_12": 0.5,
        # Chabbobboma:
        "80202_2": 2.0,
        "80202_4": 2.0,
        "80202_7": 2.0,
        "80202_8": 2.0,
        "80202_9": 2.0,
        "80202_10": 2.0,
        "80203_6": 2.0,
        "80203_7": 2.0,
        # Chisanga:
        "80204_1": 1.5,
        "80204_4": 2.0,
        "80204_6": 1.5,
        "80204_7": 1.5,
        "80204_8": 2.0,
        "80204_9": 1.5,
        "80204_10": 1.5,
        # Chiyabi:
        "81102_1": 2.0,
        "81102_2": 2.0,
        "81102_3": 2.0,
        "81102_4": 2.0,
        "81102_5": 2.0,
        "81102_6": 2.0,
        "81102_9": 2.0,
        # Luumbo:
        "80208_1": 2.0,
        "80208_4": 2.0,
        "80208_5": 2.0,
        "80208_7": 2.0,
        "80208_8": 2.0,
        "80208_10": 2.0,
        # Munyumbwe:
        "80209_2": 1.5,
        "80209_5": 0.67,
        "80209_6": 2.0,
        "80209_7": 1.5,
        "80209_8": 0.67,
        # Nyanga Chamwe:
        "80210_3": 0.15,
        "80210_4": 0.15,
        "80210_5": 0.15,
        "80210_6": 0.15,
        "80210_7": 0.15,
        "80210_8": 0.15,
        "80210_9": 0.15,
        "80210_10": 0.15,
        "80210_11": 0.15,
        "80210_12": 0.15,
        # Sinafala:
        "80204_7": 2.0,
        "80211_1": 2.0,
        "80211_3": 2.0,
        "80211_4": 2.0,
        "80211_5": 1.5,
        "80211_6": 2.0,
        "80211_7": 2.0,
        "80211_8": 1.5,
        # Sinamalima:
        "81111_1": 2.0,
        "81111_2": 2.0,
        "81111_3": 2.0,
        "81111_4": 1.5,
        "81111_5": 2.0,
        "81111_6": 2.0,
        "81111_7": 1.5,
        "81111_8": 2.0,
    }
    try:
        return [arab_param * mod_factor_dict[cl_id],
                funest_param * mod_factor_dict[cl_id]]
    except KeyError:
        return [arab_param, funest_param]

def milen_cluster_climate_category(cluster_id, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    # Return the climate category for a given cluster ID
    lookup = pd.read_csv(base + "data/milen_clusters/cluster_climate_mapping.csv")
    return lookup[lookup['cluster'] == cluster_id]['climate.cat'].item()


def closest_milen_immunity_overlay_filename(arab_param, funest_param, climate_category,
                                            base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    # Find filename of immunity file which is closest match, within given climate category, to specified arabiensis/funestus larval habitat params
    def find_nearest_from_precomputed(x, precomputed_x):
        precomp_x_arr = np.array(precomputed_x)  # assume precomputed_x is a list
        closest = np.argmin(np.abs(x - precomp_x_arr))
        return precomputed_x[closest]

    # Get arab and funest precomputed values for specified climate category
    precomputed_arab_params, precomputed_funest_params = get_arab_funest_precomputed_values_for_climate_category(
        climate_category)

    nearest_arab_param = find_nearest_from_precomputed(arab_param, precomputed_arab_params)
    nearest_funest_param = find_nearest_from_precomputed(funest_param, precomputed_funest_params)

    # Return filename for these parameters:
    fn = base + "data/immunity/{}_1_node/{}_1_node_immune_init_p1_{}_p2_{}.json".format(climate_category,
                                                                                        climate_category,
                                                                                        nearest_funest_param,
                                                                                        nearest_arab_param)
    return fn


def get_arab_funest_precomputed_values_for_climate_category(climate_category,
                                                            base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
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

    arab_params = []  # np.zeros(len(file_list))
    funest_params = []  # np.zeros(len(file_list))

    # Loop over all filenames and extract vector parameters from the filenames
    # i = 0
    # print file_list
    for fn in file_list:
        a, f = get_arab_funest_from_filename(fn)
        # arab_params[i] = a
        # funest_params[i] = f
        arab_params.append(a)
        funest_params.append(f)

    # return [np.unique(arab_params),np.unique(funest_params)]
    unique_arab_params = list(set(arab_params))
    unique_funest_params = list(set(funest_params))

    unique_arab_params.sort()
    unique_funest_params.sort()
    return [unique_arab_params, unique_funest_params]


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


def scatter_lat_long_on_map(lon, lat,
                            C=None, S=None, cbar_label=None, cmap=False, clim=None,
                            lat_range=None, lon_range=None, savefig=False, title=None,
                            cbar_or_legend='cbar'):
    # lat, lon are coordinates of data
    # C is some third data quantity used to color points
    # S is some fourth quantity used to size points
    import mpl_toolkits.basemap as base

    # fig, ax = plt.subplots()
    fig = plt.figure(figsize=(10, 10))
    ax = plt.subplot()

    if isinstance(lat_range, list) and isinstance(lon_range, list):
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

        x_min = x_min - 0.1 * w
        x_max = x_max + 0.1 * w
        y_min = y_min - 0.1 * h
        y_max = y_max + 0.1 * h

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
    if isinstance(S, int):
        S = np.ones_like(lon) * S
    elif not isinstance(S, np.ndarray):
        S = np.ones_like(lon) * 5

    if not isinstance(C, np.ndarray):
        C = 'C0'

    if not cmap:
        cmap = plt.cm.viridis

    if not isinstance(clim, list):
        clim = [min(C), max(C)]

    sc = ax.scatter(lon, lat, c=C, vmin=clim[0], vmax=clim[1], marker='s', s=S, cmap=cmap, edgecolors='black')

    if cbar_or_legend == 'cbar':
        plt.colorbar(sc)
    elif cbar_or_legend == 'legend':
        plt.legend()

    plt.title(title)

    if savefig == False:
        plt.show()
    else:
        print("saving fig to ", savefig)
        plt.tight_layout()
        plt.savefig(savefig)

    return ax


def return_satellite_map_on_plt_axes(ax, lon_range, lat_range):
    import mpl_toolkits.basemap as base

    x_min = lon_range[0]
    x_max = lon_range[1]
    y_min = lat_range[0]
    y_max = lat_range[1]

    w, h = x_max - x_min, y_max - y_min

    x_min = x_min - 0.1 * w
    x_max = x_max + 0.1 * w
    y_min = y_min - 0.1 * h
    y_max = y_max + 0.1 * h

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

    return ax


############################################################################################################

# Parse event recorder
def aggregate_events_in_recorder(recorder_df, event_type,
                                 nodeset="all",
                                 day_num_range="all",
                                 aggregate_type="week",
                                 start_date="2007-01-01"):
    # Aggregate events from the event recorder.
    # Return dictionary linking nodeID --> array of # of events, versus time

    import datetime
    from relative_time import convert_to_day_365, convert_to_date_365

    def convert_from_datestr_to_datetime(dt_str):
        dt_split = dt_str.split("-")
        return datetime.date(np.int32(dt_split[0]), np.int32(dt_split[1]), np.int32(dt_split[2]))

    def get_week_num(daynum):
        dt_string = convert_to_date_365(daynum, start_date)
        dt_datefmt = convert_from_datestr_to_datetime(dt_string)
        week_num = dt_datefmt.isocalendar()[1]
        return week_num

    def get_month_num(daynum):
        dt_string = convert_to_date_365(daynum, start_date)
        dt_datefmt = convert_from_datestr_to_datetime(dt_string)
        # print dt_string
        # print dt_datefmt
        month_num = dt_datefmt.month
        return month_num

    def get_year_num(daynum):
        dt_string = convert_to_date_365(daynum, start_date)
        dt_datefmt = convert_from_datestr_to_datetime(dt_string)
        year_num = dt_datefmt.year
        return year_num

    # print recorder_df
    event_df = recorder_df[recorder_df['Event_Name'] == event_type]
    # event_df.groupby(['Node_ID','Time']).nunique()['Time']

    events_by_day = pd.DataFrame({
        'count': event_df.groupby(["Node_ID", "Time"]).size()
    })
    events_by_day = events_by_day.reset_index()

    if nodeset != "all":
        events_by_day = events_by_day[np.in1d(events_by_day["Node_ID"], nodeset)]
    if day_num_range != "all":
        events_by_day = events_by_day[np.logical_and(events_by_day["Time"] >= day_num_range[0],
                                                     events_by_day["Time"] < day_num_range[1])]

    # print events_by_day
    # print "events by day min: ",events_by_day.min()

    if aggregate_type == "day":
        return events_by_day
    else:
        # print events_by_day
        # print "foo ",events_by_day.apply(lambda x: get_year_num(np.int32(x['Time'])), axis=1)
        events_by_day['year'] = events_by_day.apply(lambda x: get_year_num(np.int32(x['Time'])), axis=1)
        if aggregate_type == "week":
            events_by_day['week'] = events_by_day.apply(lambda x: get_week_num(np.int32(x['Time'])), axis=1)

            events_by_week = events_by_day.groupby(["Node_ID", "year", "week"]).sum()['count']
            events_by_week = events_by_week.reset_index()
            return events_by_week

        elif aggregate_type == "month":
            events_by_day['month'] = events_by_day.apply(lambda x: get_month_num(np.int32(x['Time'])), axis=1)
            events_by_month = events_by_day.groupby(["Node_ID", "year", "month"]).sum()['count']
            events_by_month = events_by_month.reset_index()
            return events_by_month


############################################################################################################

# HFCA_additional_CHW_lookup = {
#     "Chabbobboma": ["so Chipepo RHC Siancheka",
#                     "so Chipepo RHC Chilindi",
#                     "so Gulumunyanga Health Post (CHA)",
#                     "so Gulumunyanga HP Hamatuba"],
#
#     "Munyumbwe": ["so Makuyu HP Ganikoongo",
#                   "so Makuyu HP Katete",
#                   "so Fumbo  Delivery Point (step 2b)"],
#
# }


HFCA_milen_cluster_lookup = {
    "Bbondo": ["80201_1",
               "80201_2",
               "80201_3",
               "80201_4",
               "80201_5",
               "80201_6",
               "80201_7",
               "80201_8",
               "80201_9",
               "80201_10",
               "80201_11",
               "80201_12"],
    "Chabbobboma": ["80202_1",
                    "80202_2",
                    "80202_3",
                    "80202_4",
                    "80202_5",
                    "80202_6",
                    "80202_7",
                    "80202_8",
                    "80202_9",
                    "80202_10",
                    "80203_4",
                    "80203_6",
                    "80203_7",
                    "80208_7"],
    "Chisanga": ["80204_1",
                 "80204_2",
                 "80204_3",
                 "80204_4",
                 "80204_5",
                 "80204_6",
                 "80204_7",
                 "80204_8",
                 "80204_9",
                 "80204_10"],
    "Chiyabi": ["81102_1",
                "81102_2",
                "81102_3",
                "81102_4",
                "81102_5",
                "81102_6",
                "81102_7",
                "81102_8",
                "81102_9"],
    "Luumbo": ["80202_1",
               "80208_1",
               "80208_2",
               "80208_3",
               "80208_4",
               "80208_5",
               "80208_6",
               "80208_7",
               "80208_8",
               "80208_9",
               "80208_10",
               "80210_3"],
    "Munyumbwe": ["80209_2",
                  "80209_3",
                  "80209_4",
                  "80209_5",
                  "80209_6",
                  "80209_7",
                  "80209_8",
                  "80209_9",
                  "80209_10"],
    "Nyanga Chaamwe": ["80210_1",
                       "80210_3",
                       "80210_4",
                       "80210_5",
                       "80210_6",
                       "80210_7",
                       "80210_8",
                       "80210_9",
                       "80210_10",
                       "80210_11",
                       "80210_12"],
    "Sinafala": ["80204_7",
                 "80211_1",
                 "80211_2",
                 "80211_3",
                 "80211_4",
                 "80211_5",
                 "80211_6",
                 "80211_7",
                 "80211_8"],
    "Sinamalima": ["81111_1",
                   "81111_2",
                   "81111_3",
                   "81111_4",
                   "81111_5",
                   "81111_6",
                   "81111_7",
                   "81111_8"]
}



############################################################################################################

def compute_round_date(round,
                       cell_ids,
                       weight="covpop",
                       start_date="2007-01-01",
                       base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/',
                       prev_fn="data/prevalence/2017-12-20/raw/grid_prevalence_with_dates.csv",
                       pop_fn="data/gridded_pop/cleaned/all_max_pop.csv"):
    # Open Caitlin's file
    prev_df = pd.read_csv(base + prev_fn)
    prev_df = prev_df[prev_df['round']==round]

    df = pd.DataFrame({
        "cell_ids": cell_ids
    })

    df = df.merge(prev_df,how='left',left_on="cell_ids",right_on="grid_cell")
    # Drop NANs (which occur when a given cell doesn't appear in this round)
    df = df.dropna()

    if len(df) == 0:
        return -1
    df['sim_day'] = df.apply(lambda x: convert_to_day_365(x['date'],start_date),axis=1)

    # Remove possible outlier dates
    # non_outlier_dates = nonoutlier_mask(df['sim_day'])
    # if np.sum(non_outlier_dates) < 2:
    #     return -1
    # else:
    #     df = df[non_outlier_dates]

    # Merge in population information, so that we can weight properly.
    if weight == "covpop":
        # pop_df = pd.read_csv(base + "data/interventions/kariba/2017-11-27/raw/grid_prevalence.csv")
        pop_df = pd.read_csv(base + prev_fn)
        pop_df = pop_df[pop_df["round"]==round]

        df = df.merge(pop_df,how="left",left_on="cell_ids",right_on="grid_cell")
        weighted_round_date = int((df["N_x"]*df["sim_day"]).sum()/df["N_x"].sum())


    elif weight == "maxpop":
        pop_df = pd.read_csv(base + pop_fn)

        df = df.merge(pop_df, how="left", left_on="cell_ids", right_on="node_label")
        weighted_round_date = int((df["pop"] * df["sim_day"]).sum() / df["pop"].sum())

    print("Weighted round date: ",convert_to_date_365(weighted_round_date,start_date))
    return weighted_round_date

def round_date_sanity_check():
    for rd in range(1,11):
        compute_round_date(rd,np.arange(10000),weight="covpop")
        compute_round_date(rd,np.arange(10000),weight="maxpop")

def nonoutlier_mask(data, m=2.):
    # https://stackoverflow.com/questions/11686720/is-there-a-numpy-builtin-to-reject-outliers-from-a-list
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d / mdev if mdev else 0.
    return s < m



############################################################################################################

def add_cell_intervention_timing_rugs_to_plot(ax,
                                              cell_ids,
                                              start_date="2007-01-01",
                                              base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/',
                                              irs_relative_path="data/interventions/kariba/2017-11-27/raw/grid_all_irs_events.csv",
                                              itn_relative_path="data/interventions/kariba/2017-11-27/raw/grid_all_itn_events.csv",
                                              mda_relative_path="data/interventions/kariba/2017-11-27/raw/grid_all_mda_events.csv",
                                              ymax=1.0):
    import matplotlib.dates as mdates
    import seaborn as sns
    sns.set_style("darkgrid")

    # start_date = "2007-01-01"  # Day 1 of simulation
    date_format = "%Y-%m-%d"

    foo = mdates.strpdate2num(date_format)


    # Plot vertical lines for different intervention timepoints:
    # IRS:
    irs_df = pd.read_csv(base + irs_relative_path)
    irs_df = irs_df[np.in1d(np.array(irs_df["grid_cell"]), cell_ids)]
    print("plotting IRS lines")
    lbl_flag = 0
    for d in irs_df['fulldate']:
        if lbl_flag == 0:
            lbl = "IRS events"
            lbl_flag = 1
        else:
            lbl = None
        ax.axvline(foo(d), c='C0', ymin=ymax*0.8, ymax=ymax, lw=0.5, alpha=0.4, label=lbl, zorder=1)
    print("done plotting IRS lines")

    # ITN:
    itn_df = pd.read_csv(base + itn_relative_path)
    itn_df = itn_df[np.in1d(np.array(itn_df["grid_cell"]), cell_ids)]
    print("plotting itn lines")
    lbl_flag = 0
    for d in itn_df['fulldate']:
        if lbl_flag == 0:
            lbl = "ITN events"
            lbl_flag = 1
        else:
            lbl = None
        ax.axvline(foo(d), c='C1', ymin=ymax*0.5, ymax=ymax*0.7, lw=0.5, alpha=0.4, label=lbl, zorder=1)
    print("done plotting itn lines")

    # MDA:
    mda_df = pd.read_csv(base + mda_relative_path)
    mda_df = mda_df[np.in1d(np.array(mda_df["grid_cell"]), cell_ids)]
    print("plotting mda lines")
    lbl_flag = 0
    for d in mda_df['fulldate']:
        if lbl_flag == 0:
            lbl = "MDA events"
            lbl_flag = 1
        else:
            lbl = None
        ax.axvline(foo(d), c='C2', ymin=0, ymax=ymax*0.2, lw=0.5, alpha=0.2, label=lbl, zorder=2)
    print("done plotting mda lines")

    # # MSAT:
    # msat_df = pd.read_csv(base + "data/interventions/kariba/2017-11-27/raw/grid_all_msat_events.csv")
    # msat_df = msat_df[np.in1d(np.array(msat_df["grid_cell"]), cell_ids)]
    # print ("plotting msat lines")
    # lbl_flag = 0
    # for d in msat_df['fulldate']:
    #     if lbl_flag == 0:
    #         lbl = "MSAT events"
    #         lbl_flag = 1
    #     else:
    #         lbl = None
    #     ax.axvline(foo(d), c='C3', ymin=0, ymax=0.2, lw=0.5, alpha=0.2, label=lbl, zorder=2)
    # print("done plotting msat lines")

    # STEPD:
    # stepd_df = pd.read_csv(base + "data/interventions/kariba/2017-11-27/raw/grid_all_stepd_events.csv")
    # stepd_df = stepd_df[np.in1d(np.array(stepd_df["grid_cell"]), cell_ids)]
    # print ("plotting stepd lines")
    # lbl_flag = 0
    # for d in stepd_df['fulldate']:
    #     if lbl_flag == 0:
    #         lbl = "CHWs added"
    #         lbl_flag = 1
    #     else:
    #         lbl = None
    #     ax.axvline(foo(d), c='C4', ymin=0, ymax=1.0, lw=1.0, label=lbl, linestyle='dashed', zorder=2)
    # print("done plotting stepd lines")



def find_cells_for_this_milen_cluster(mc_name, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
    mc_lookup_df = pd.read_csv(base + "data/milen_clusters/cluster_to_grid_lookup.csv")
    in_mc = mc_lookup_df['cluster_id'] == mc_name
    cell_list = list(set(mc_lookup_df['grid_cell'][in_mc]))
    cell_list.sort()
    return np.array(cell_list)



############################################################################################################
# Calibration/Analyzer helper functions

# def load_



def safe_start_year_duration_for_climate_generator(start_year,sim_duration,fixed_year=2014):
    # Return a "safe" start year for the climate generator.
    # Conditions: start year must be >= 2001, and need to make sure that the climate for the fixed_year is correct

    cg_duration = np.min([sim_duration,10])


    if start_year > 2001:
        return [start_year,cg_duration]
    elif (fixed_year-start_year) <= cg_duration:
        return [start_year, cg_duration]
    else:
        for i in np.flip(np.arange(11),axis=0):
            if (fixed_year-start_year) % i == 0:
                return [fixed_year-i,i]





############################################################################################################
# Mozambique calibration:

