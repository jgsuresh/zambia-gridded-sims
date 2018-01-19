from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from relative_time import *
from simtools.AnalyzeManager.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
import numpy as np
import pandas as pd

class RDTPrevAnalyzer(BaseAnalyzer):

    filenames = ['output\SpatialReport_Population.bin', 'output\SpatialReport_Prevalence.bin', 'output\SpatialReport_New_Diagnostic_Prevalence.bin']

    def __init__(self):
        super(RDTPrevAnalyzer, self).__init__()
        self.my_data = {}
        self.metadata = {}

        self.RDT_prev_by_node = {}
        self.pop_by_node = {}
        self.catch = {}

    def filter(self, sim_metadata):
        return sim_metadata['Run_Number'] == 0

    def apply(self, parser):
        exp_name = parser.experiment.exp_name
        self.catch[parser.sim_id] = exp_name.split('_')[0] # Assumes the experiment name is "CATCHNAME_full"

        pop_data = parser.raw_data[self.filenames[0]]
        prev_data = parser.raw_data[self.filenames[1]]
        RDT_prev_data = parser.raw_data[self.filenames[2]]

        self.n_tstep = pop_data['n_tstep']
        self.n_nodes = pop_data['n_nodes']


        self.RDT_prev_by_node[parser.sim_id] = {}
        self.pop_by_node[parser.sim_id] = {}
        # Initialize node arrays:
        for j in range(self.n_nodes):
            self.RDT_prev_by_node[parser.sim_id][j] = np.zeros(self.n_tstep)
            self.pop_by_node[parser.sim_id][j] = np.zeros(self.n_tstep)

        # Collect per-node data:
        for i in range(self.n_tstep):
            RDT_timestep_data = RDT_prev_data['data'][i]
            pop_timestep_data = pop_data['data'][i]
            for j in range(self.n_nodes):
                self.RDT_prev_by_node[parser.sim_id][j][i] = RDT_timestep_data[j]
                self.pop_by_node[parser.sim_id][j][i] = pop_timestep_data[j]

    def finalize(self):
        print("")

    def plot(self):
        import matplotlib.pyplot as plt
        from matplotlib import cm
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

        print(daydates_mdates)

        maxpop_ever = 0
        for sim_id, data in self.RDT_prev_by_node.items():
            for j in range(self.n_nodes):
                maxp = np.max(self.pop_by_node[sim_id][j])
                if maxp > maxpop_ever:
                    maxpop_ever = maxp

        plt.figure(figsize=(12,5))
        for sim_id, data in self.RDT_prev_by_node.items():
            for j in range(self.n_nodes):
                maxp = np.max(self.pop_by_node[sim_id][j])
                pop_ratio = np.float(maxp)/np.float(maxpop_ever)
                plt.plot_date(daydates_mdates, self.RDT_prev_by_node[sim_id][j],fmt='-',alpha=pop_ratio,color=cm.Greens_r(pop_ratio),zorder=1)

        catch = self.catch.itervalues().next()

        # Look up catchment prevalence data from precomputed file:
        base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
        df = pd.read_csv(base + "data/interventions/kariba/2017-11-27/cleaned/catch_prevalence.csv")
        catch_prev = np.array(df[catch])

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
            plt.scatter(round_dates_array[:-4], catch_prev[:-4], c='red', s=70, label='Data for {}'.format(catch.capitalize()),zorder=10)
            plt.scatter(round_dates_array[-4:], catch_prev[-4:], c='gray', s=70 ,label='HFCA not in MDA round',zorder=10)
        elif catch in ["chisanga"]:
            plt.scatter(np.append(round_dates_array[:3],round_dates_array[5:]), np.append(catch_prev[:3],catch_prev[5:]), c='red', s=70, label='Data for {}'.format(catch.capitalize()),zorder=10)
            plt.scatter(round_dates_array[3:5], catch_prev[3:5], c='gray', s=70, label='HFCA not in MDA round',zorder=10)
        else:
            plt.scatter(round_dates_array, catch_prev, c='red', s=70, label='Data for {}'.format(catch.capitalize()),zorder=10)






        # # Plot Chiyabi interventions as vertical lines:
        # if True:
        #     # plot_only_first = False
        #
        #     # Event information files
        #     itn_event_file = "C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/cbever/chiyabi_hfca_itn_events.csv"
        #     irs_event_file = "C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/cbever/chiyabi_hfca_irs_events.csv"
        #     msat_event_file = "C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/cbever/chiyabi_hfca_msat_events.csv"
        #     mda_event_file = "C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/cbever/chiyabi_hfca_mda_events.csv"
        #     healthseek_event_file = "C:/Users/jsuresh//OneDrive - IDMOD/Code/zambia/cbever/chiyabi_hfca_healthseek_events.csv"
        #     stepd_event_file = "C:/Users/jsuresh//OneDrive - IDMOD/Code/zambia/cbever/chiyabi_hfca_stepd_events.csv"
        #
        #     # Import event info
        #     itn_events = pd.read_csv(itn_event_file)
        #     irs_events = pd.read_csv(irs_event_file)
        #     msat_events = pd.read_csv(msat_event_file)
        #     mda_events = pd.read_csv(mda_event_file)
        #     healthseek_events = pd.read_csv(healthseek_event_file)
        #     stepd_events = pd.read_csv(stepd_event_file)
        #
        #     # for date in itn_events['fulldate']: plt.axvline(convert_to_day(date, start_date, date_format=date_format), ls='dashed', color='C0') #,label='ITN Events')
        #     # for date in irs_events['fulldate']: plt.axvline(convert_to_day(date, start_date, date_format=date_format), ls='dashed', color='C1') #,label='IRS Events')
        #     # for date in msat_events['fulldate']: plt.axvline(convert_to_day(date, start_date, date_format=date_format), ls='dashed', color='C2') #,label='MSAT Events')
        #     # for date in mda_events['fulldate']: plt.axvline(convert_to_day(date, start_date, date_format=date_format), ls='dashed', color='C3') #,label='MDA Events')
        #
        #     for date in itn_events['fulldate']: plt.axvline(foo(date), ls='dashed', color='blue') #,label='ITN Events')
        #     for date in irs_events['fulldate']: plt.axvline(foo(date), ls='dashed', color='green') #,label='IRS Events')
        #     for date in msat_events['fulldate']:plt.axvline(foo(date), ls='dashed', color='red') #,label='MSAT Events')
        #     for date in mda_events['fulldate']: plt.axvline(foo(date), ls='dashed', color='purple') #,label='MDA Events')
        #
        #
        #     # colors:
        #     # c0 = blue = ITN
        #     # c1 = green = IRS
        #     # c2 = red = MSAT
        #     # c3 = purple = MDA




        plt.legend()
        # plt.xlim([3000,7000])
        plt.xlim([foo("2010-01-01"), foo("2019-01-01")])
        # plt.show()
        plt.tight_layout()
        plt.savefig(base + "data/figs/{}_prev_node.png".format(catch))


if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()

    # am.add_experiment(retrieve_experiment("43cac760-cbd6-e711-9414-f0921c16b9e5")) # bbondo
    # am.add_experiment(retrieve_experiment("a31b516a-cbd6-e711-9414-f0921c16b9e5"))  # chabbobboma
    # am.add_experiment(retrieve_experiment("1ecdf372-cbd6-e711-9414-f0921c16b9e5")) # chisanga
    # am.add_experiment(retrieve_experiment("957e6159-32d6-e711-9414-f0921c16b9e5")) # chiyabi
    # am.add_experiment(retrieve_experiment("9669907b-cbd6-e711-9414-f0921c16b9e5"))  # luumbo
    # am.add_experiment(retrieve_experiment("fbe40809-ccd6-e711-9414-f0921c16b9e5"))  # munyumbwe
    am.add_experiment(retrieve_experiment("8aadd6a0-cbd6-e711-9414-f0921c16b9e5"))  # nyanga chaamwe
    # am.add_experiment(retrieve_experiment("d18a9aa8-cbd6-e711-9414-f0921c16b9e5"))  # sinafala
    # am.add_experiment(retrieve_experiment("d28a9aa8-cbd6-e711-9414-f0921c16b9e5"))  # sinamalima


    am.add_analyzer(RDTPrevAnalyzer())
    am.analyze()