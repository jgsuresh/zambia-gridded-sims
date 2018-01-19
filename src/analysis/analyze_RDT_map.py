from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from relative_time import *
from simtools.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
import numpy as np
import pandas as pd

from gridded_sim_general import *
from RDT_map_from_data import *
from datetime import date

class RDTPrevAnalyzer(BaseAnalyzer):

    filenames = ['output/SpatialReport_Population.bin', 'output/SpatialReport_New_Diagnostic_Prevalence.bin', 'Assets/Demographics/demo.json']

    def __init__(self):
        super(RDTPrevAnalyzer, self).__init__()
        self.my_data = {}
        self.metadata = {}

        self.RDT_prev_by_node = {}
        self.pop_by_node = {}
        self.catch = {}
        self.demo_dict = None
        self.node_ids = {}

        self.base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'

    def filter(self, sim_metadata):
        # return True
        return sim_metadata['Run_Number'] == 0

    def apply(self, parser):
        exp_name = parser.experiment.exp_name
        self.catch = exp_name.split('_')[0] # Assumes the experiment name is "CATCHNAME_full"

        pop_data = parser.raw_data[self.filenames[0]]
        RDT_prev_data = parser.raw_data[self.filenames[1]]
        self.demo_dict = parser.raw_data[self.filenames[2]]

        self.node_ids = pop_data['nodeids']
        self.n_tstep = pop_data['n_tstep']
        self.n_nodes = pop_data['n_nodes']

        self.RDT_prev_by_node = {}
        self.pop_by_node = {}
        # Initialize node arrays:
        for j in range(self.n_nodes):
            self.RDT_prev_by_node[j] = np.zeros(self.n_tstep)
            self.pop_by_node[j] = np.zeros(self.n_tstep)

        # Collect per-node data:
        for i in range(self.n_tstep):
            RDT_timestep_data = RDT_prev_data['data'][i]
            pop_timestep_data = pop_data['data'][i]
            for j in range(self.n_nodes):
                self.RDT_prev_by_node[j][i] = RDT_timestep_data[j]
                self.pop_by_node[j][i] = pop_timestep_data[j]

    def finalize(self):
        pass


    def plot(self):
        import matplotlib.pyplot as plt
        from matplotlib import cm
        import matplotlib.dates as mdates
        from matplotlib import ticker
        import seaborn as sns
        sns.set_style("darkgrid")

        def convert_date_str_to_date_obj(date_str):
            # Break into pieces:
            date_list = date_str.split("-")
            for i in range(3):
                date_list[i] = int(date_list[i])
            return date(date_list[0],date_list[1],date_list[2])

        start_date = "2007-01-01"  # Day 1 of simulation
        date_format = "%Y-%m-%d"
        cmap = plt.get_cmap("RdBu",11)
        mode = "ratio" #"diff"

        foo = mdates.strpdate2num(date_format)

        # Convert simulation day numbers to actual dates:
        daynum = np.arange(self.n_tstep)
        daydates_list = []
        daydates_mdates = np.array([])
        for dayn in daynum:
            hold = convert_to_date(dayn, start_date, date_format=date_format)
            daydates_list.append(hold)
            daydates_mdates = np.append(daydates_mdates,foo(hold))

        # Look up catchment prevalence data from precomputed file:
        obs_df = pd.read_csv(self.base + "data/interventions/kariba/2017-11-27/raw/grid_prevalence.csv")

        #fixme Get round dates from prevalence data.
        round_dates = ["2012-07-01","2012-09-30","2012-11-30","2013-07-01","2013-08-31","2013-10-31","2014-12-31","2015-03-01","2015-09-30","2016-02-29"]


        # Get get grid cell IDs, and lat/long for each DTK node:
        node_grid_cell_ids = convert_from_dtk_node_ids_to_grid_cells_using_demo(self.node_ids, self.demo_dict)
        [node_lat, node_lon] = get_lat_long_dtk_nodes(self.node_ids, self.demo_dict)

        # sim_df = pd.DataFrame({
        #     "node_lat": node_lat,
        #     "node_lon": node_lon,
        #     "sim_grid_cell_ids": node_grid_cell_ids
        # })

        # Define custom colorbar limits to match already-generated figures:
        clim_dict = {
            'bbondo': 0.2,
            'chabbobboma': 0.6,
            'chisanga': 0.5,
            'chiyabi': 0.7,
            'luumbo': 0.8,
            'munyumbwe': 0.5,
            'nyanga chaamwe': 0.2,
            'sinafala': 0.6,
            'sinamalima': 0.7
        }

        # Loop over all rounds:
        rd = 1
        for rd_date in round_dates:
            print("Working on round {}".format(rd))

            # Get the simulation-day-number for this particular round
            sd = convert_date_str_to_date_obj(start_date)
            ed = convert_date_str_to_date_obj(rd_date)
            day_num = (ed-sd).days

            # Find the node data for this date:
            pop = np.zeros(self.n_nodes)
            prev = np.zeros(self.n_nodes)

            for j in range(self.n_nodes):
                pop[j] = self.pop_by_node[j][day_num]
                prev[j] = self.RDT_prev_by_node[j][day_num]

            # Create dataframe of simulation data in prep for merging with observational data:
            data_df = pd.DataFrame({
                "grid_cell_id": node_grid_cell_ids,
                "sim_pop": pop,
                "sim_prev": prev,
            })

            # Get the observational data for this round:
            obs_round_df = obs_df[obs_df['round']==rd]

            # Merge the two:
            data_df = data_df.merge(obs_round_df,how='left',left_on='grid_cell_id',right_on='grid_cell')
            data_df = data_df.rename(columns={'prev': 'obs_prev', 'N': 'obs_pop'})
            data_df['obs_pop'] = data_df['obs_pop'].fillna(0)


            # Plot:
            plt.figure(figsize=(10, 10))
            ax = plt.subplot(111)

            ax = return_satellite_map_on_plt_axes(ax,
                                                  [np.min(node_lon), np.max(node_lon)],
                                                  [np.min(node_lat), np.max(node_lat)])

            S = scale_pt_size(pop, size_min=50, size_max=500)

            # Plot bad points, all in gray:
            missing_from_obs = data_df['obs_pop'] < (1./5.) * data_df['sim_pop']
            if np.sum(missing_from_obs) > 0:
                ax.scatter(node_lon[missing_from_obs],
                           node_lat[missing_from_obs],
                           s=S[missing_from_obs],
                           c='gray',
                           marker='s',
                           edgecolor='black',
                           label='Obs data missing or unreliable')

            # Plot good points, colored by difference in prevalence:
            in_obs = np.logical_not(missing_from_obs)
            if np.sum(in_obs) > 0:
                if mode == "diff":
                    prev_diff = (data_df[in_obs]['sim_prev']-data_df[in_obs]['obs_prev'])#/(data_df[in_obs]['obs_prev'])


                scale = max([np.max(prev),np.abs(np.min(prev))])

                sc = ax.scatter(node_lon[in_obs],
                                node_lat[in_obs],
                                s=S[in_obs],
                                c=prev_diff,
                                cmap=cmap,
                                vmin=-clim_dict[self.catch],
                                vmax=clim_dict[self.catch],
                                marker='s',
                                edgecolor='black',
                                # vmin=clim[0], vmax=clim[1],
                                )

                cb=plt.colorbar(sc)
                tick_locator = ticker.MaxNLocator(nbins=10)
                cb.locator = tick_locator
                cb.update_ticks()
                cb.set_label('(Sim Prev) - (Measured Prev)')

            ax.set_title('SIM: {} round {}'.format(self.catch.capitalize(), rd))
            # plt.show()
            plt.savefig(self.base + "data/figs/RDT_diff_maps/{}_rd_{}_diff.png".format(self.catch, rd))

            plt.close('all')
            rd += 1





        # round_dates_mdate = []
        # for i in range(10):
        #     day_mdate = foo(round_dates[i])
        #     round_dates_mdate.append(day_mdate)
        # round_dates_array = np.array(round_dates_mdate)
        #
        # [node_lat, node_lon] = get_lat_long_dtk_nodes(node_ids)
        #
        # rd = 1
        # for rd_date in round_dates_array:


        # scatter_lat_long_on_map(node_lon,node_lat,C=)


        # plt.legend()
        # # plt.xlim([3000,7000])
        # plt.xlim([foo("2010-01-01"), foo("2019-01-01")])
        # plt.show()
        # plt.tight_layout()
        # plt.savefig(base + "data/figs/{}_prev_node.png".format(catch))


if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager.AnalyzeManager()

    # am.add_experiment(retrieve_experiment("43cac760-cbd6-e711-9414-f0921c16b9e5")) # bbondo
    # am.add_experiment(retrieve_experiment("a31b516a-cbd6-e711-9414-f0921c16b9e5"))  # chabbobboma
    # am.add_experiment(retrieve_experiment("1ecdf372-cbd6-e711-9414-f0921c16b9e5")) # chisanga
    # am.add_experiment(retrieve_experiment("957e6159-32d6-e711-9414-f0921c16b9e5")) # chiyabi
    # am.add_experiment(retrieve_experiment("9669907b-cbd6-e711-9414-f0921c16b9e5"))  # luumbo
    # am.add_experiment(retrieve_experiment("fbe40809-ccd6-e711-9414-f0921c16b9e5"))  # munyumbwe
    # am.add_experiment(retrieve_experiment("8aadd6a0-cbd6-e711-9414-f0921c16b9e5"))  # nyanga chaamwe
    # am.add_experiment(retrieve_experiment("d18a9aa8-cbd6-e711-9414-f0921c16b9e5"))  # sinafala
    am.add_experiment(retrieve_experiment("d28a9aa8-cbd6-e711-9414-f0921c16b9e5"))  # sinamalima


    am.add_analyzer(RDTPrevAnalyzer())
    am.analyze()