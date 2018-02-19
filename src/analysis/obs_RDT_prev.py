import pandas as pd
import numpy as np

# Get list of all catchments
# For each catchment, find all grid cells which correspond to this catchment
# Do a pop-weighted sum of the prevalence in this grid cell, at each available round
# Save the results in a new intermediate file



# Clean up grid lookup file by making sure the catchment names are capitalized:
# def write_new_lookup_with_caps_names(grid_lookup):
#     # Outdated as of 11/27/17: Caitlin's new script generates all catchment names in lower case:
#     # Fix problem where only some of the catchment names in the lookup file are capitalized:
#     def capitalize_if_possible(catch):
#         if type(catch) == float and np.isnan(catch):
#             return catch
#         else:
#             return catch.capitalize()
#
#     catch_arr = np.array(grid_lookup['catchment'])
#     catch_arr_cap = map(capitalize_if_possible,catch_arr)
#
#     grid_lookup_clean = grid_lookup.copy()
#     grid_lookup_clean['catchment'] = pd.Series(catch_arr_cap,index=grid_lookup_clean.index)
#     # grid_prev['catchment'] = grid_lookup.applymap(lambda x: capitalize_if_possible(x['catchment']))
#
#     grid_lookup_clean.to_csv(base + 'data/interventions/kariba/2017-11-21/cleaned/grid_lookup_caps.csv')

def get_catch_list():
    # grid_lookup = pd.read_csv(base + 'data/interventions/kariba/2017-11-21/cleaned/grid_lookup_caps.csv')
    # grid_lookup = pd.read_csv(base + 'data/interventions/kariba/2017-11-27/raw/grid_lookup.csv')
    grid_lookup = pd.read_csv(base + 'data/mozambique/cleaned/grid_lookup.csv')
    # Get a list of all catchments:
    return grid_lookup['catchment'].unique()

def get_grid_cells_for_catch(grid_lookup,catch):
    cells = grid_lookup['grid_cell'][grid_lookup['catchment']==catch]
    return np.array(cells)

def get_catch_prev(full_df,catch,weighting='population'):
    catch_cells = full_df['grid_cell'][full_df['catchment']==catch]
    in_catch = np.in1d(grid_prev['grid_cell'],np.array(catch_cells))

    catch_prev = {}
    denom = {}
    for rd in range(1,11):
        in_rd = full_df['round'] == rd
        full_cut = np.logical_and(in_catch,in_rd)
        if np.sum(full_cut) > 0:
            in_catch_this_rd = full_df[np.logical_and(in_catch,in_rd)]

            prev = np.array(in_catch_this_rd['prev'])
            pop_this_rd = np.array(in_catch_this_rd['N'])
            maxpop_ever = np.array(in_catch_this_rd['pop'])

            # Exclude nan pixels:
            good_pix = np.logical_not(np.isnan(prev))
            prev = prev[good_pix]
            pop_this_rd = pop_this_rd[good_pix]
            maxpop_ever = maxpop_ever[good_pix]

            if np.sum(good_pix) > 0:

                if weighting=='population':
                    total_prev = np.float(np.sum(prev*maxpop_ever))/np.float(np.sum(maxpop_ever))
                    catch_prev[rd] = total_prev
                    denom[rd] = np.sum(maxpop_ever)

                elif weighting=='MDA_coverage':
                    total_prev = np.float(np.sum(prev*pop_this_rd))/np.float(np.sum(pop_this_rd))
                    catch_prev[rd] = total_prev
                    denom[rd] = np.sum(pop_this_rd)

    return [catch_prev,denom]

def save_prev_by_catch(full_df, weighting='population'):
    # For each catchment, find all grid cells which correspond to this catchment
    catch_list = get_catch_list()
    prev_dict = {}
    denom_dict = {}
    for catch in catch_list:
        # Ignore the NAN - pixels which are not in a catchment:
        if type(catch)== str and catch != "NONE":
            [prev_dict[catch],denom_dict[catch]] = get_catch_prev(full_df,catch,weighting=weighting)

    # Turn this into a dataframe, and save it to CSV
    prev_df = pd.DataFrame(prev_dict)
    denom_df = pd.DataFrame(denom_dict)
    # prev_df = pd.concat([prev_df,denom_df],axis=1)
    # prev_df = pd.DataFrame({
    #     "prev": prev_dict,
    #     "weight": denom_dict
    # })

    if weighting == 'population':
        # prev_df.to_csv(base + 'data/interventions/kariba/2017-11-27/cleaned/catch_prevalence_pop_weighted.csv')
        prev_df.to_csv(base + 'data/mozambique/cleaned/catch_prevalence_pop_weighted.csv')
        denom_df.to_csv(base + 'data/mozambique/cleaned/catch_prevalence_pop_weighted_denom.csv')
    elif weighting == 'MDA_coverage':
        # prev_df.to_csv(base + 'data/interventions/kariba/2017-11-27/cleaned/catch_prevalence_coverage_weighted.csv')
        prev_df.to_csv(base + 'data/mozambique/cleaned/catch_prevalence_coverage_weighted.csv')
        denom_df.to_csv(base + 'data/mozambique/cleaned/catch_prevalence_coverage_weighted_denom.csv')



if __name__ == "__main__":
    base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
    # grid_prev = pd.read_csv(base + 'data/interventions/kariba/2017-11-27/raw/grid_prevalence.csv')
    # grid_lookup = pd.read_csv(base + 'data/interventions/kariba/2017-11-27/raw/grid_lookup.csv')
    grid_prev = pd.read_csv(base + 'data/mozambique/grid_prevalence_with_dates.csv')
    grid_lookup = pd.read_csv(base + 'data/mozambique/cleaned/grid_lookup.csv')

    # Merge so each grid cell also knows its own catchment
    full_df = grid_prev.merge(grid_lookup,how='left',left_on='grid_cell',right_on='grid_cell')

    # Merge again so each grid cell also knows its own maximum population
    grid_maxpop = pd.read_csv(base + 'data/mozambique/grid_population.csv')
    full_df = full_df.merge(grid_maxpop,how='left',left_on='grid_cell',right_on='node_label')

    # write_new_lookup_with_caps_names(grid_prev)
    # grid_lookup = pd.read_csv(base + 'data/interventions/kariba/2017-11-27/cleaned/grid_lookup_caps.csv')

    save_prev_by_catch(full_df, weighting='population')
    save_prev_by_catch(full_df, weighting='MDA_coverage')

