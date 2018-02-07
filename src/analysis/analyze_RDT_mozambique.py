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

        self.base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'

    def filter(self, sim_metadata):
        # if sim_metadata["sim_id"] == "08dff6bf-690b-e811-9415-f0921c16b9e5":
        return True
        # else:
        #     return False

    def apply(self, parser):
        exp_name = parser.experiment.exp_name
        # self.catch[parser.sim_id] = exp_name.split('_')[0] # Assumes the experiment name is "CATCHNAME_full"
        self.catch[parser.sim_id] = "Motaze" #"Mahel"

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
        print("")

    def plot(self):
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates
        import seaborn as sns
        sns.set_style("darkgrid")

        # start_date = "2007-01-01"  # Day 1 of simulation
        start_date = "1994-01-01"
        date_format = "%Y-%m-%d"

        foo = mdates.strpdate2num(date_format)

        daynum = np.arange(self.n_tstep)
        daydates_list = []
        daydates_mdates = np.array([])
        for dayn in daynum:
            hold = convert_to_date(dayn, start_date, date_format=date_format)
            daydates_list.append(hold)
            daydates_mdates = np.append(daydates_mdates,foo(hold))

        # print(daydates_mdates)

        plt.figure(figsize=(12,5))
        ax = plt.subplot(111)

        lbl_flag = 0
        for sim_id, data in self.RDT_prev_aggr.items():
            if lbl_flag == 0:
                lbl = "Simulations"
                lbl_flag = 1
            else:
                lbl = None
            plt.plot_date(daydates_mdates, self.RDT_prev_aggr[sim_id],fmt='-',color='black',label=lbl,lw=1.2,zorder=10)

        # catch = self.catch.values().next()
        catch = list(self.catch.values())[0]
        # catch = catch.lower()

        catch_cell_ids = find_cells_for_this_catchment(catch,path_from_base='data/mozambique/grid_lookup.csv')

        ###############################################################################################################

        ###############################################################################################################

        # Look up catchment prevalence data from precomputed file:
        df = pd.read_csv(self.base + "data/mozambique/cleaned/catch_prevalence_coverage_weighted.csv")
        catch_prev_cov_weighted = np.array(df[catch])
        df = pd.read_csv(self.base + "data/mozambique/cleaned/catch_prevalence_pop_weighted.csv")
        catch_prev_pop_weighted = np.array(df[catch])


        # global_round_dates = ["2012-06-18", "2012-08-29", "2012-11-03", "2013-06-09", "2013-08-11", "2013-10-08",
        #                        "2014-12-17", "2015-02-17", "2015-09-20", "2016-02-16"]
        global_round_dates = ["2000-01-01"]*5
        # global_round_dates is computed from gridded_sim_general.round_date_sanity_check()

        # Compute round dates that are appropriate for the given catchment

        # rd_day_sim = np.zeros(10)
        round_dates = list(global_round_dates)
        for rd in range(1,6):
            rd_day_sim = compute_round_date(rd,catch_cell_ids,
                                            start_date=start_date,
                                            weight="maxpop",
                                            prev_fn="data/mozambique/grid_prevalence_with_dates.csv",
                                            pop_fn="data/mozambique/grid_population.csv"
                                            )

            if rd_day_sim != -1:
                # Then computing round date failed.  Default to global round date in this case:
                round_dates[rd-1] = convert_to_date_365(rd_day_sim, start_date)
            # rd_day_sim[rd-1] =



        round_dates_mdate = []
        for i in range(5):
            day_mdate = foo(round_dates[i])
            round_dates_mdate.append(day_mdate)
        round_dates_array = np.array(round_dates_mdate)

        # if catch in ["chabbobboma","chipepo","gwembe","lukande","nyanga chaamwe"]:
        #     plt.scatter(round_dates_array[:-4], catch_prev_cov_weighted[:-4], c='pink',edgecolors='black', s=70, alpha=0.9,label='{}: Coverage-weighted RDT+'.format(catch.capitalize()),zorder=20)
        #     plt.scatter(round_dates_array[-4:], catch_prev_cov_weighted[-4:], c='gray',edgecolors='black', s=70,zorder=20) # ,label='HFCA not in MDA round'
        # elif catch in ["chisanga"]:
        #     plt.scatter(np.append(round_dates_array[:3],round_dates_array[5:]), np.append(catch_prev_cov_weighted[:3],catch_prev_cov_weighted[5:]), c='pink',edgecolors='black', s=70, alpha=0.9,label='{}: Coverage-weighted RDT+'.format(catch.capitalize()),zorder=20)
        #     plt.scatter(round_dates_array[3:5], catch_prev_cov_weighted[3:5], c='gray',edgecolors='black', s=70,zorder=20) #, label='HFCA not in MDA round'
        # else:
        #     plt.scatter(round_dates_array, catch_prev_cov_weighted, c='pink',edgecolors='black', s=70, alpha=0.9,label='{}: Coverage-weighted RDT+'.format(catch.capitalize()),zorder=20)
        #
        # if catch in ["chabbobboma","chipepo","gwembe","lukande","nyanga chaamwe"]:
        #     plt.scatter(round_dates_array[:-4], catch_prev_pop_weighted[:-4], c='teal',edgecolors='black', marker='s',s=70, alpha=0.8,label='{}: Pop-weighted RDT+'.format(catch.capitalize()),zorder=20)
        #     plt.scatter(round_dates_array[-4:], catch_prev_pop_weighted[-4:], c='gray',edgecolors='black', marker='s',s=70,zorder=20) # ,label='HFCA not in MDA round'
        # elif catch in ["chisanga"]:
        #     plt.scatter(np.append(round_dates_array[:3],round_dates_array[5:]), np.append(catch_prev_pop_weighted[:3],catch_prev_pop_weighted[5:]), c='teal',edgecolors='black', marker='s',s=70, alpha=0.8,label='{}: Pop-weighted RDT+'.format(catch.capitalize()),zorder=20)
        #     plt.scatter(round_dates_array[3:5], catch_prev_pop_weighted[3:5], c='gray',edgecolors='black', marker='s',s=70,zorder=20) # , label='HFCA not in MDA round'
        # else:
        plt.scatter(round_dates_array, catch_prev_pop_weighted, c='teal',edgecolors='black', marker='s', s=70, alpha=0.8,label='{}: Pop-weighted RDT+'.format(catch.capitalize()),zorder=20)

        # For each round time-point, also plot the corresponding MDA-coverage-weighted (not full pop-weighted) RDT prevalence from the simulation:

        # Get the observational data for how many people were observed in that round, in each pixel.
        # Divide this by the "max pop ever seen in this pixel" to get the "MDA coverage" for that pixel.
        # Aggregate an MDA-coverage-weighted RDT prevalence from the corresponding pixels in the simulation.
        prev_df = pd.read_csv(self.base + "data/mozambique/grid_prevalence_with_dates.csv")
        max_pop_df = pd.read_csv(self.base + "data/mozambique/grid_population.csv")

        full_df = prev_df.merge(max_pop_df,how='left',left_on='grid_cell',right_on='node_label')


        # catch_df = full_df[np.in1d(full_df['grid_cell'],catch_cell_ids)]

        # # Loop over every round
        # coverage_corrected_prev_sim = np.zeros(10)
        # for round in range(1,10):
        #     in_round = catch_df['round'] == round
        #     temp_df = catch_df[in_round]
        #
        #     # For each round, get list of cells, their populations, and their MDA-coverage's
        #     # Find way to convert from cell ID list to node ID list.
        #     # Find way to get node ID list from



        add_cell_intervention_timing_rugs_to_plot(ax,
                                                  catch_cell_ids,
                                                  start_date=start_date,
                                                  irs_relative_path="data/mozambique/grid_all_irs_events.csv",
                                                  itn_relative_path="data/mozambique/grid_all_itn_events.csv",
                                                  mda_relative_path="data/mozambique/grid_all_mda_events.csv",
                                                  ymax=1.0)

        plt.legend()
        # plt.xlim([3000,7000])
        plt.xlim([foo("2015-01-01"), foo("2018-01-01")])
        plt.ylim([-0.01,0.1])

        plt.tight_layout()
        plt.show()
        # plt.savefig(self.base + "data/figs/{}_prev.png".format(catch))


if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()

    # Calibration experiments:
    # am.add_experiment(retrieve_experiment("09829129-b00b-e811-9415-f0921c16b9e5")) #Mahel
    am.add_experiment(retrieve_experiment("8f552695-c90b-e811-9415-f0921c16b9e5")) #Motaze

    am.add_analyzer(RDTPrevAnalyzer())
    am.analyze()