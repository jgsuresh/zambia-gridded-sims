from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from relative_time import *
from simtools.AnalyzeManager.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
import numpy as np
import pandas as pd
from gridded_sim_general import *

class RDTPrevAnalyzer(BaseAnalyzer):

    filenames = ['output/SpatialReport_Population.bin', 'output/SpatialReport_Prevalence.bin', 'output/SpatialReport_New_Diagnostic_Prevalence.bin', 'Assets/Demographics/demo.json']

    def __init__(self):
        super(RDTPrevAnalyzer, self).__init__()
        self.my_data = {}
        self.metadata = {}

        self.prev_aggr = {}
        self.RDT_prev_aggr = {}
        self.catch = {}
        self.node_ids = {}

    def filter(self, sim_metadata):
        return True

    def apply(self, parser):
        exp_name = parser.experiment.exp_name
        self.catch[parser.sim_id] = exp_name.split('_')[0] # Assumes the experiment name is "CATCHNAME_full"

        pop_data = parser.raw_data[self.filenames[0]]
        prev_data = parser.raw_data[self.filenames[1]]
        RDT_prev_data = parser.raw_data[self.filenames[2]]
        demo = parser.raw_data[self.filenames[3]]

        self.node_ids[parser.sim_id] = pop_data['nodeids']
        self.n_tstep = pop_data['n_tstep']
        self.n_nodes = pop_data['n_nodes']

        # Get initial population of nodes:
        self.pop_init = np.zeros(self.n_nodes)
        for ni in range(self.n_nodes):
            self.pop_init[ni] = pop_data['data'][0][ni]

        # Collect aggregated data:
        self.prev_aggr[parser.sim_id] = np.zeros(self.n_tstep)
        self.RDT_prev_aggr[parser.sim_id] = np.zeros(self.n_tstep)
        for i in range(self.n_tstep):
            self.prev_aggr[parser.sim_id][i] = np.sum(pop_data['data'][i]*prev_data['data'][i])/np.sum(pop_data['data'][i])
            self.RDT_prev_aggr[parser.sim_id][i] = np.sum(pop_data['data'][i] * RDT_prev_data['data'][i]) / np.sum(pop_data['data'][i])

    def finalize(self):
        print ""

    def plot(self):
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        import seaborn as sns
        sns.set_style("darkgrid")

        start_date = "2007-01-01"  # Day 1 of simulation
        date_format = "%Y-%m-%d"

        foo = mdates.strpdate2num(date_format)

        daynum = np.arange(self.n_tstep)
        daydates_list = []
        daydates_mdates = np.array([])
        for dayn in daynum:
            hold = convert_to_date(dayn, start_date, date_format=date_format)
            daydates_list.append(hold)
            daydates_mdates = np.append(daydates_mdates,foo(hold))

        print daydates_mdates

        plt.figure(figsize=(12,5))
        for sim_id, data in self.RDT_prev_aggr.items():
            plt.plot_date(daydates_mdates, self.RDT_prev_aggr[sim_id],fmt='-')

        catch = self.catch.itervalues().next()

        # Look up catchment prevalence data from precomputed file:
        base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
        df = pd.read_csv(base + "data/interventions/kariba/2017-11-27/cleaned/catch_prevalence_coverage_weighted.csv")
        catch_prev_cov_weighted = np.array(df[catch])
        df = pd.read_csv(base + "data/interventions/kariba/2017-11-27/cleaned/catch_prevalence_pop_weighted.csv")
        catch_prev_pop_weighted = np.array(df[catch])


        round_dates = ["2012-07-01","2012-09-30","2012-11-30","2013-07-01","2013-08-31","2013-10-31","2014-12-31","2015-03-01","2015-09-30","2016-02-29"]
        # round_dates = {
        #     "1": "2012-07-01",
        #     "2": "2012-09-30",
        #     "3": "2012-11-30",
        #     "4": "2013-07-01",
        #     "5": "2013-08-31",
        #     "6": "2013-10-31",
        #     "7": "2014-12-31",
        #     "8": "2015-03-01",
        #     "9": "2015-09-30",
        #     "10":"2016-02-29"
        # }

        round_dates_mdate = []
        for i in range(10):
            day_mdate = foo(round_dates[i])
            round_dates_mdate.append(day_mdate)
        round_dates_array = np.array(round_dates_mdate)

        if catch in ["chabbobboma","chipepo","gwembe","lukande","nyanga chaamwe"]:
            plt.scatter(round_dates_array[:-4], catch_prev_cov_weighted[:-4], c='red', s=70, alpha=0.5,label='{}: Coverage-weighted RDT+'.format(catch.capitalize()))
            plt.scatter(round_dates_array[-4:], catch_prev_cov_weighted[-4:], c='gray', s=70) # ,label='HFCA not in MDA round'
        elif catch in ["chisanga"]:
            plt.scatter(np.append(round_dates_array[:3],round_dates_array[5:]), np.append(catch_prev_cov_weighted[:3],catch_prev_cov_weighted[5:]), c='red', s=70, alpha=0.5,label='{}: Coverage-weighted RDT+'.format(catch.capitalize()))
            plt.scatter(round_dates_array[3:5], catch_prev_cov_weighted[3:5], c='gray', s=70) #, label='HFCA not in MDA round'
        else:
            plt.scatter(round_dates_array, catch_prev_cov_weighted, c='red', s=70, alpha=0.5,label='{}: Coverage-weighted RDT+'.format(catch.capitalize()))

        if catch in ["chabbobboma","chipepo","gwembe","lukande","nyanga chaamwe"]:
            plt.scatter(round_dates_array[:-4], catch_prev_pop_weighted[:-4], c='blue', marker='s',s=70, alpha=0.5,label='{}: Pop-weighted RDT+'.format(catch.capitalize()))
            plt.scatter(round_dates_array[-4:], catch_prev_pop_weighted[-4:], c='gray', marker='s',s=70) # ,label='HFCA not in MDA round'
        elif catch in ["chisanga"]:
            plt.scatter(np.append(round_dates_array[:3],round_dates_array[5:]), np.append(catch_prev_pop_weighted[:3],catch_prev_pop_weighted[5:]), c='blue', marker='s',s=70, alpha=0.5,label='{}: Pop-weighted RDT+'.format(catch.capitalize()))
            plt.scatter(round_dates_array[3:5], catch_prev_pop_weighted[3:5], c='gray', marker='s',s=70) # , label='HFCA not in MDA round'
        else:
            plt.scatter(round_dates_array, catch_prev_pop_weighted, c='blue', marker='s', s=70, alpha=0.5,label='{}: Pop-weighted RDT+'.format(catch.capitalize()))

        # For each round time-point, also plot the corresponding MDA-coverage-weighted (not full pop-weighted) RDT prevalence from the simulation:

        # Get the observational data for how many people were observed in that round, in each pixel.
        # Divide this by the "max pop ever seen in this pixel" to get the "MDA coverage" for that pixel.
        # Aggregate an MDA-coverage-weighted RDT prevalence from the corresponding pixels in the simulation.
        prev_df = pd.read_csv(base + "data/interventions/kariba/2017-11-27/raw/grid_prevalence.csv")
        max_pop_df = pd.read_csv(base + "data/gridded_pop/cleaned/all_max_pop.csv")

        full_df = prev_df.merge(max_pop_df,how='left',left_on='grid_cell',right_on='node_label')

        catch_cell_ids = find_cells_for_this_catchment(catch)
        catch_df = full_df[np.in1d(full_df['grid_cell'],catch_cell_ids)]

        # # Loop over every round
        # coverage_corrected_prev_sim = np.zeros(10)
        # for round in range(1,10):
        #     in_round = catch_df['round'] == round
        #     temp_df = catch_df[in_round]
        #
        #     # For each round, get list of cells, their populations, and their MDA-coverage's
        #     # Find way to convert from cell ID list to node ID list.
        #     # Find way to get node ID list from




        plt.legend()
        # plt.xlim([3000,7000])
        plt.xlim([foo("2010-01-01"), foo("2019-01-01")])
        # plt.show()
        plt.tight_layout()
        plt.savefig(base + "data/figs/{}_prev.png".format(catch))


if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()


    # Corrected vector param immunity, BUT incorrect StepD (only ~1 cell got stepd)
    # am.add_experiment(retrieve_experiment("5539241e-32d6-e711-9414-f0921c16b9e5")) # bbondo
    # am.add_experiment(retrieve_experiment("15f5282b-32d6-e711-9414-f0921c16b9e5"))  # chabbobboma
    # am.add_experiment(retrieve_experiment("80169448-32d6-e711-9414-f0921c16b9e5")) # chisanga
    # am.add_experiment(retrieve_experiment("957e6159-32d6-e711-9414-f0921c16b9e5")) # chiyabi
    # am.add_experiment(retrieve_experiment("7aa7e969-32d6-e711-9414-f0921c16b9e5"))  # luumbo
    # am.add_experiment(retrieve_experiment("695f9d80-32d6-e711-9414-f0921c16b9e5"))  # munyumbwe
    # am.add_experiment(retrieve_experiment("1ea21996-32d6-e711-9414-f0921c16b9e5"))  # nyanga chaamwe
    # am.add_experiment(retrieve_experiment("eff77db9-32d6-e711-9414-f0921c16b9e5"))  # sinafala
    # am.add_experiment(retrieve_experiment("451fe5d4-32d6-e711-9414-f0921c16b9e5"))  # sinamalima

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

    am.add_analyzer(RDTPrevAnalyzer())
    am.analyze()