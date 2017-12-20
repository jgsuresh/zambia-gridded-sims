import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")

from gridded_sim_general import *

base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
df = pd.read_csv(base + "data/incidence/CHW_grid_lookup.csv")

catch_list = ['bbondo','chabbobboma','chisanga','chiyabi','luumbo','munyumbwe','nyanga chaamwe','sinafala','sinamalima']

for catch in catch_list[6:]:
#     catch = catch_list[5]
    grid_cells = find_cells_for_this_catchment(catch)

    catch_df = pd.DataFrame({
        'grid_cell_ids': grid_cells
    })

    catch_df = catch_df.merge(df,how='left',left_on='grid_cell_ids',right_on='loc.id')


    # Add "chw_flag" for plotting purposes
    chw_names = list(set(catch_df['closest.chw.name']))
    num_chw = len(chw_names)
    chw_ind_df = pd.DataFrame({
        "chw_name": chw_names,
        "chw_flag": range(num_chw)
    })

    catch_df = catch_df.merge(chw_ind_df,how='left',left_on='closest.chw.name',right_on='chw_name')

    plt.figure(figsize=(10,10))
    ax = plt.subplot(111)

    ax = return_satellite_map_on_plt_axes(ax,
                                         [np.min(catch_df['mid.x']), np.max(catch_df['mid.x'])],
                                         [np.min(catch_df['mid.y']), np.max(catch_df['mid.y'])])

    chw_color = 0
    outsider_color = 0
    m_list = ['o','s','^']
    for chw_name in chw_names:
        this_chw = catch_df['closest.chw.name'] == chw_name

        if catch.title() not in chw_name:
            marker='x'
            cnum = outsider_color
            outsider_color += 1
        else:
            marker = m_list[chw_color % 3]
            cnum = chw_color
            chw_color += 1

        ax.scatter(catch_df[this_chw]['mid.x'],catch_df[this_chw]['mid.y'],label=chw_name,
                   c='C{}'.format(cnum % 8),
                   marker=marker,
                   s=150)



    # Scatter plot the positions of the actual CHWs themselves:
    chw_locs_df = pd.read_csv(base + "data/incidence/MACEPA ZM DHIS2 Org Units Only Date Format Fix CHW placement v2.csv")
    chw_ind_df = chw_ind_df.merge(chw_locs_df,how='left',left_on='chw_name',right_on='Org.Unit.Name')
    ax.scatter(chw_ind_df['Longitude'],chw_ind_df['Latitude'],label="CHW Locations",
                   c='yellow'.format(chw_color % 8),
                   edgecolor='black',
                   marker='*',
                   s=200)


    # Scatter plot positions of health facilities:
    #

    ax.set_title(catch.title())

    legend=ax.legend(frameon=True,framealpha=0.5,facecolor='white',fontsize=8)
    # legend.get_frame().set_facecolor('white')
    plt.show()
    # scatter_lat_long_on_map(catch_df['mid.x'],catch_df['mid.y'],C=catch_df['chw_flag'])