import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")

from gridded_sim_general import *


def scale_pt_size(x,size_min=2,size_max=100):
    if np.max(x) != np.min(x):
        a = np.float(size_max-size_min)/np.float(np.max(x)-np.min(x))
        b = size_min - a*np.min(x)
        b = np.min([b,0])
        return a*x+b
    else:
        return np.ones_like(x) * np.average([size_min, size_max])

def cbar_scale(prev):
    import math
    def roundup(x, v=10.0):
        return int(math.ceil(x / v)) * int(v)

    prev = np.nan_to_num(prev)
    cmax = np.percentile(prev,80)
    if cmax > 0.1:
        cmax = roundup(cmax*100.)/100.
    elif cmax < 0.1:
        cmax = roundup(cmax*100.,v=5.0)/100.
    return [0,cmax]


if __name__ == "__main__":
    # Make a 2D map of RDT prevalence in a catchment:
    # - One map for each round, with title showing round # and date
    # - Color of pixel corresponds to RDT prevalence (include colorbar)
    # - Size of pixel corresponds to population
    # - Background shows geography
    # Import data:
    base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
    prev_df = pd.read_csv(base + "data/interventions/kariba/2017-11-27/raw/grid_prevalence.csv")
    lookup_df = pd.read_csv(base + "data/interventions/kariba/2017-11-27/raw/grid_lookup.csv")

    # Merge the two databases so that we can directly get the lat and long for each cell:
    full_df = prev_df.merge(lookup_df,how='left',left_on='grid_cell',right_on='grid_cell')

    # Plot-specific stuff:
    cmap = plt.get_cmap('Blues',5)


    # milen_catch_list = ['bbondo','chabbobboma','chisanga','chiyabi','luumbo','munyumbwe','nyanga chaamwe','sinafala','sinamalima']
    milen_catch_list = ['chabbobboma']

    for catch in milen_catch_list:
        cell_ids = find_cells_for_this_catchment(catch)
        in_catch = np.in1d(full_df['grid_cell'], cell_ids)

        lon_range = [np.min(full_df[in_catch]['mid_x']), np.max(full_df[in_catch]['mid_x'])]
        lat_range = [np.min(full_df[in_catch]['mid_y']), np.max(full_df[in_catch]['mid_y'])]

        w, h = lon_range[1]-lon_range[0], lat_range[1]-lat_range[0]

        lon_range[0] = lon_range[0] - 0.1 * w
        lon_range[1] = lon_range[1] + 0.1 * w
        lat_range[0] = lat_range[0] - 0.1 * h
        lat_range[1] = lat_range[1] + 0.1 * h


        cflag = 0
        for round in range(1,11):
        # round = 10
            in_rd = full_df['round'] == round

            in_catch_rd = np.logical_and(in_catch,in_rd)

            if np.sum(in_catch_rd) > 1:
                temp_df = full_df[in_catch_rd]
                S = scale_pt_size(np.array(temp_df['N']),size_min=50,size_max=500)
                C = np.array(temp_df['prev'])

                if cflag == 0:
                    clim = cbar_scale(C)
                    cflag = 1

                fname = base + "data/figs/RDT_maps/wiki/{}_rd_{}_data.png".format(catch,round)

                scatter_lat_long_on_map(np.array(temp_df['mid_x']),
                                        np.array(temp_df['mid_y']),
                                        C=np.array(temp_df['prev']),
                                        S=S,
                                        lat_range=lat_range,
                                        lon_range=lon_range,
                                        clim=clim,
                                        cbar_label='Prevalence',
                                        cmap=cmap,
                                        title='DATA: {} round {}'.format(catch.capitalize(),round),
                                        savefig=fname)

                plt.close('all')