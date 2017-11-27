# From survey data, create gridding and make a max-population-over-all-rounds

import pandas as pd
import numpy as np


def get_catch_grid(catch_name):
    base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
    grid_df = pd.read_csv(base + "data/gridded_pop/raw/grid_lookup.csv")
    if catch_name == 'all':
        catch_df = grid_df
    else:
        print "catch_name = ",catch_name
        catch_df = grid_df[grid_df['catchment']==catch_name]
    catch_df = catch_df.reset_index()
    return catch_df

def get_latlong_lists_grid(catch_df):
    def find_min_nonzero_diff(x):
        x = np.sort(x)
        dx = np.abs(x[1:]-x[:-1])
        min_nonzero_dx = np.min(dx[dx > 0])
        return min_nonzero_dx

    dx = find_min_nonzero_diff(np.array(catch_df['mid_x']))
    dy = find_min_nonzero_diff(np.array(catch_df['mid_y']))

    x_min = np.array(catch_df['mid_x']) - dx/2.
    x_max = np.array(catch_df['mid_x']) + dx/2.
    y_min = np.array(catch_df['mid_y']) - dy/2.
    y_max = np.array(catch_df['mid_y']) + dy/2.
    return [x_min, x_max, y_min, y_max]

def find_max_pix_pop(survey_df, lat_bnds, long_bnds):
    lat_cut = np.logical_and(survey_df['latitude'] > lat_bnds[0],
                             survey_df['latitude'] < lat_bnds[1])
    long_cut = np.logical_and(survey_df['longitude'] > long_bnds[0],
                              survey_df['longitude'] < long_bnds[1])
    in_pix = np.logical_and(lat_cut,long_cut)
    pix_df = survey_df[in_pix]

    rd_pop = np.zeros(10)
    for i in xrange(10):
        rd_pop[i] = np.sum(pix_df['round']==i+1)

    return np.max(rd_pop)

def save_in_old_format(catch_df,outname):
    # Rename other columns to match previous format:
    catch_df = catch_df.drop('index',axis=1)
    catch_df = catch_df.drop('catchment',axis=1)
    catch_df = catch_df.rename(columns={'grid_cell':'node_label','mid_x': 'lon', 'mid_y': 'lat'})

    catch_df.to_csv(outname)
    return catch_df


# Loop over each grid cell, and find max pop over all rounds:
def compute_max_pop_catch_cells(catch_name):
    base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
    survey_df = pd.read_csv(base + "data/raw/masterDatasetAllRounds2012-2016.csv")

    catch_df = get_catch_grid(catch_name)
    long_min,long_max,lat_min,lat_max = get_latlong_lists_grid(catch_df)

    n_pix = len(long_min)
    max_pop = np.zeros(n_pix)
    for i in np.arange(n_pix):
        print "On pixel {} of {}...".format(i+1,n_pix)
        lat_bnds = [lat_min[i],lat_max[i]]
        long_bnds = [long_min[i],long_max[i]]
        max_pop[i] = find_max_pix_pop(survey_df,lat_bnds,long_bnds)

    catch_df['pop'] = pd.Series(max_pop,index=catch_df.index)

    outname = base + "data/gridded_pop/cleaned/{}_max_pop.csv".format(catch_name.lower())
    catch_df = save_in_old_format(catch_df,outname)
    return catch_df

def gen_pop_csv_files_for_milen_catchments():
    # milen_catch_list = ['chisanga','sinamalima','chiyabi','sinafala','munyumbwe','chipepo','chipepo siavonga','lukonde',
    #                     'chabbobboma','bbondo','luumbo','nyanga chaamwe']
    milen_catch_list = ['chisanga','sinamalima','chiyabi','sinafala','munyumbwe','chipepo','chipepo siavonga',
                        'chabbobboma','bbondo','luumbo','nyanga chaamwe']


    for catch in milen_catch_list:
        compute_max_pop_catch_cells(catch)


if __name__== "__main__":
    gen_pop_csv_files_for_milen_catchments()
# cdf = compute_max_pop_catch_cells("all")
# cdf = compute_max_pop_catch_cells("Mapatizya")
# cdf = compute_max_pop_catch_cells("Munyumbwe")


