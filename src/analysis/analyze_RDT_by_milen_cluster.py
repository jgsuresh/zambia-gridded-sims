from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from relative_time import *
from simtools.AnalyzeManager.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
import numpy as np
import pandas as pd
from gridded_sim_general import *

class RDTPrevAnalyzer(BaseAnalyzer):

    filenames = ['output/SpatialReport_Population.bin', 'output/SpatialReport_New_Diagnostic_Prevalence.bin', 'Assets/Demographics/demo.json']

    def __init__(self):
        super(RDTPrevAnalyzer, self).__init__()
        self.my_data = {}
        self.metadata = {}

        self.RDT_prev_by_node = {}
        self.catch = {}
        self.node_ids = {}

        self.base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'

    def filter(self, sim_metadata):
        # return True
        return sim_metadata['Run_Number'] == 0

    def apply(self, parser):
        exp_name = parser.experiment.exp_name
        self.catch = exp_name.split('_')[0] # Assumes the experiment name is "CATCHNAME_other"

        pop_data = parser.raw_data[self.filenames[0]]
        RDT_prev_data = parser.raw_data[self.filenames[1]]
        self.demo = parser.raw_data[self.filenames[2]]

        self.node_ids = pop_data['nodeids']
        self.n_tstep = pop_data['n_tstep']
        self.n_nodes = pop_data['n_nodes']

        # Get initial population of nodes:
        self.pop_init = np.zeros(self.n_nodes)
        for ni in range(self.n_nodes):
            self.pop_init[ni] = pop_data['data'][0][ni]

        # Collect aggregated data:
        self.RDT_prev_by_node = {}

        for j in range(self.n_nodes):
            node_id = self.node_ids[j]
            self.RDT_prev_by_node[node_id] = np.zeros(self.n_tstep)

            for i in range(self.n_tstep):
                self.RDT_prev_by_node[node_id][i] = RDT_prev_data['data'][i][j]
                # self.RDT_prev_aggr[parser.sim_id][i] = np.sum(pop_data['data'][i] * RDT_prev_data['data'][i]) / np.sum(pop_data['data'][i])

    def finalize(self):
        print ""

    def plot(self):
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        import seaborn as sns
        sns.set_style("darkgrid")

        start_date = "2007-01-01"  # Day 1 of simulation
        date_format = "%Y-%m-%d"

        # Get cells for this catchment
        catch_cell_ids = find_cells_for_this_catchment(self.catch)

        # Get actual dates for the simulation days:
        foo = mdates.strpdate2num(date_format)

        daynum = np.arange(self.n_tstep)
        daydates_list = []
        daydates_mdates = np.array([])
        for dayn in daynum:
            hold = convert_to_date(dayn, start_date, date_format=date_format)
            daydates_list.append(hold)
            daydates_mdates = np.append(daydates_mdates,foo(hold))

        print daydates_mdates



        plt.figure(figsize=(12, 5))

        # Get lookup files ready:
        mc_lookup_df = pd.read_csv(self.base + "data/milen_clusters/cluster_lookup.csv") #fixme
        prev_lookup_df = pd.read_csv() #fixme

        # Loop over every Milen-cluster that is associated with this HFCA (from dictionary in gridded_sim_general.py)
        num_mc = len(HFCA_milen_cluster_lookup[self.catch])

        for mi in range(num_mc):
            mc = HFCA_milen_cluster_lookup[self.catch][mi]

            ax = plt.subplot(num_mc / 2, 2, mi + 1)
            ax.set_title(mc)

            # For each Milen-cluster, find the cells which correspond to this Milen-cluster (from one of Caitlin's lookup files)
            cells_this_mc = mc_lookup_df['grid'][] #fixme

            # Convert these cells into node IDs:
            node_ids = convert_from_grid_cells_to_dtk_node_ids_using_demo(cells_this_mc, self.demo)

            # For each of these nodes, plot:
            #   - the simulated RDT prevalence of that particular node
            #   - the observed RDT prevalence points (from one of Caitlin's lookup files), sized according to the population of the node.
            for ni in range(len(node_ids)):
                node_id = node_ids[ni]
                grid_cell = cells_this_mc[ni]
                pop = self.pop_init[ni] #fixme

                # Plot simulated RDT prevalence for this node:
                ax.plot_date(daydates_mdates,self.RDT_prev_by_node[node_id]
                             ,c='C{}'.format(8%ni))

                # Plot observed RDT prevalence for corresponding grid cell:
                this_node = prev_lookup_df['grid'] == grid_cell
                ax.plot_date(prev_lookup_df['date'][this_node],prev_lookup_df['prev'][this_node],linestyle='none',
                             c = 'C{}'.format(8 % ni),s=np.sqrt(pop))

            ax.set_xlabel("Date")
            ax.set_ylabel("RDT Prevalence")
            ax.legend()
            ax.set_xlim([foo("2010-01-01"), foo("2019-01-01")])

        plt.tight_layout()
        plt.show()
        # plt.savefig(self.base + "data/figs/{}_prev_by_mcluster.png".format(catch))


if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()

    # Corrected stepd
    # am.add_experiment(retrieve_experiment("43cac760-cbd6-e711-9414-f0921c16b9e5")) # bbondo
    # am.add_experiment(retrieve_experiment("a31b516a-cbd6-e711-9414-f0921c16b9e5"))  # chabbobboma
    # am.add_experiment(retrieve_experiment("1ecdf372-cbd6-e711-9414-f0921c16b9e5")) # chisanga
    # am.add_experiment(retrieve_experiment("957e6159-32d6-e711-9414-f0921c16b9e5")) # chiyabi
    # am.add_experiment(retrieve_experiment("9669907b-cbd6-e711-9414-f0921c16b9e5"))  # luumbo
    # am.add_experiment(retrieve_experiment("fbe40809-ccd6-e711-9414-f0921c16b9e5"))  # munyumbwe
    # am.add_experiment(retrieve_experiment("8aadd6a0-cbd6-e711-9414-f0921c16b9e5"))  # nyanga chaamwe
    # am.add_experiment(retrieve_experiment("d18a9aa8-cbd6-e711-9414-f0921c16b9e5"))  # sinafala
    am.add_experiment(retrieve_experiment("d28a9aa8-cbd6-e711-9414-f0921c16b9e5"))  # sinamalima

    # Old MBGSR
    # am.add_experiment(retrieve_experiment("7f188957-2fe1-e711-9414-f0921c16b9e5")) # bbondo
    # am.add_experiment(retrieve_experiment("f60d69eb-2fe1-e711-9414-f0921c16b9e5"))  # chabbobboma
    # am.add_experiment(retrieve_experiment("7aa30068-2fe1-e711-9414-f0921c16b9e5")) # chisanga
    # am.add_experiment(retrieve_experiment("d57bccae-25e1-e711-9414-f0921c16b9e5")) # chiyabi
    # am.add_experiment(retrieve_experiment("5d5cff6d-2fe1-e711-9414-f0921c16b9e5"))  # luumbo
    # am.add_experiment(retrieve_experiment("cf37cd7b-2fe1-e711-9414-f0921c16b9e5"))  # munyumbwe
    # am.add_experiment(retrieve_experiment("94aa85fb-2fe1-e711-9414-f0921c16b9e5"))  # nyanga chaamwe
    # am.add_experiment(retrieve_experiment("f5c0fb13-30e1-e711-9414-f0921c16b9e5"))  # sinafala
    # am.add_experiment(retrieve_experiment("33b92b39-30e1-e711-9414-f0921c16b9e5"))  # sinamalima

    am.add_analyzer(RDTPrevAnalyzer())
    am.analyze()