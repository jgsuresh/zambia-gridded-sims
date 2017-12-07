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

    # filenames = ['output\SpatialReport_Population.bin',
    #              'output\SpatialReport_Prevalence.bin',
    #              'output\SpatialReport_New_Diagnostic_Prevalence.bin',
    #              'Assets\Demographics\demo.json']
    filenames = ['output/SpatialReport_Population.bin', 'output/SpatialReport_Prevalence.bin',
                 'output/SpatialReport_New_Diagnostic_Prevalence.bin', 'Assets/Demographics/demo.json']

    def __init__(self):
        super(RDTPrevAnalyzer, self).__init__()
        self.my_data = {}
        self.metadata = {}

        self.RDT_prev_by_node = {}
        self.pop_by_node = {}
        self.catch = {}
        self.demo_file = {}
        self.node_ids = {}

    def filter(self, sim_metadata):
        return sim_metadata['Run_Number'] == 0

    def apply(self, parser):
        exp_name = parser.experiment.exp_name
        self.catch = exp_name.split('_')[0] # Assumes the experiment name is "CATCHNAME_full"

        base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
        self.demo_file = base + "data/COMPS_experiments/{}_full/Demographics/demo.json".format(self.catch)

        pop_data = parser.raw_data[self.filenames[0]]
        prev_data = parser.raw_data[self.filenames[1]]
        RDT_prev_data = parser.raw_data[self.filenames[2]]
        # demo_data = parser.raw_data[self.filenames[3]]

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
        print ""

    def plot(self):
        import matplotlib.pyplot as plt
        from matplotlib import cm
        import matplotlib.dates as mdates
        import seaborn as sns
        sns.set_style("darkgrid")

        start_date = "2007-01-01"  # Day 1 of simulation
        date_format = "%Y-%m-%d"

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
        base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
        df = pd.read_csv(base + "data/interventions/kariba/2017-11-27/cleaned/catch_prevalence.csv")
        catch_prev = np.array(df[self.catch])

        round_dates = ["2012-07-01","2012-09-30","2012-11-30","2013-07-01","2013-08-31","2013-10-31","2014-12-31","2015-03-01","2015-09-30","2016-02-29"]

        def convert_date_str_to_date_obj(date_str):
            # Break into pieces:
            date_list = date_str.split("-")
            for i in range(3):
                date_list[i] = int(date_list[i])
            return date(date_list[0],date_list[1],date_list[2])

        # def get_all_node_RDT_prev_for_given_daynum(daynum):
        #
        #     self.RDT_prev_by_node

        [node_lat, node_lon] = get_lat_long_dtk_nodes(self.node_ids, self.demo_file)

        cmap = plt.get_cmap('Blues', 5)

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


        rd = 1
        for rd_date in round_dates:
            sd = convert_date_str_to_date_obj(start_date)
            ed = convert_date_str_to_date_obj(rd_date)

            day_num = (ed-sd).days

            # Find the node data for this date:
            pop = np.zeros(self.n_nodes)
            prev = np.zeros(self.n_nodes)

            for j in range(self.n_nodes):
                pop[j] = self.pop_by_node[j][day_num]
                prev[j] = self.RDT_prev_by_node[j][day_num]

            S = scale_pt_size(pop, size_min=50, size_max=500)
            C = prev
            # clim = cbar_scale(C)
            clim = [0.,clim_dict[self.catch]]
            fname = base + "data/figs/RDT_maps/sims/{}_rd_{}.png".format(self.catch, rd)

            scatter_lat_long_on_map(node_lon,
                                    node_lat,
                                    C=C,
                                    S=S,
                                    clim=clim,
                                    cbar_label='Prevalence',
                                    cmap=cmap,
                                    title='SIM: {} round {}'.format(self.catch.capitalize(), rd),
                                    savefig=fname)
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
        # # plt.show()
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
    am.add_experiment(retrieve_experiment("fbe40809-ccd6-e711-9414-f0921c16b9e5"))  # munyumbwe
    # am.add_experiment(retrieve_experiment("8aadd6a0-cbd6-e711-9414-f0921c16b9e5"))  # nyanga chaamwe
    # am.add_experiment(retrieve_experiment("d18a9aa8-cbd6-e711-9414-f0921c16b9e5"))  # sinafala
    # am.add_experiment(retrieve_experiment("d28a9aa8-cbd6-e711-9414-f0921c16b9e5"))  # sinamalima


    am.add_analyzer(RDTPrevAnalyzer())
    am.analyze()