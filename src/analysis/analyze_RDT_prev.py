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

        self.prev_aggr = {}
        self.RDT_prev_aggr = {}
        self.RDT_prev_by_node = {}

        # self.plot_labels = ["Multinode x_Local_Migration=10","Multinode x_Local_Migration=1","Multinode x_Local_Migration=0.1","Multinode x_Local_Migration=0", "Singlenode"]

    def filter(self, sim_metadata):
        # return sim_metadata['Run_Number'] == 0
        # print sim_metadata.keys()
        # print sim_metadata['run_number']
        # return sim_metadata['x_Local_Migration']==0.1
        # return parser.sim_data['Run_Number'] == 0
        # return True
        return (sim_metadata['IRS'] and sim_metadata['ITNs'] and sim_metadata['MDA'] and sim_metadata['MSAT'] and sim_metadata['StepD'] and sim_metadata['Healthseek'])

    def apply(self, parser):
        pop_data = parser.raw_data[self.filenames[0]]
        prev_data = parser.raw_data[self.filenames[1]]
        RDT_prev_data = parser.raw_data[self.filenames[2]]

        # c = parser.raw_data[self.filenames[0]].keys()

        self.n_tstep = pop_data['n_tstep']
        self.n_nodes = pop_data['n_nodes']

        # Get initial population of nodes:
        self.pop_init = np.zeros(self.n_nodes)
        for ni in range(self.n_nodes):
            self.pop_init[ni] = pop_data['data'][0][ni]
        # print self.pop_init

        # self.my_data[parser.sim_id] = data['data'][500]
        foo = pop_data['data'][:]
        foo2 = pop_data['data']
        prev_data_aggr = (pop_data['data']*prev_data['data'])#/self

        # Collect aggregated data:
        self.prev_aggr[parser.sim_id] = np.zeros(self.n_tstep)
        self.RDT_prev_aggr[parser.sim_id] = np.zeros(self.n_tstep)
        for i in range(self.n_tstep):
            self.prev_aggr[parser.sim_id][i] = np.sum(pop_data['data'][i]*prev_data['data'][i])/np.sum(pop_data['data'][i])
            self.RDT_prev_aggr[parser.sim_id][i] = np.sum(pop_data['data'][i] * RDT_prev_data['data'][i]) / np.sum(pop_data['data'][i])

        # Collect node-by-node data:
        # self.RDT_prev_by_node[parser.sim_id] = {}
        # if self.n_nodes > 1:
        #     for ni in range(self.n_nodes):
        #         self.RDT_prev_by_node[parser.sim_id][ni] = np.zeros(self.n_tstep)
        #         for i in range(self.n_tstep):
        #             self.RDT_prev_by_node[parser.sim_id][ni][i] = RDT_prev_data['data'][i][ni]
        #             # self.RDT_prev_aggr[parser.sim_id][i] = np.sum(pop_data['data'][i] * RDT_prev_data['data'][i]) / np.sum(pop_data['data'][i])
        #
        # else:
        #     self.RDT_prev_by_node[parser.sim_id]




        # self.metadata[parser.sim_id] = parser.sim_data
        # print parser.simulation.experiment.exp_name
        # print parser.sim_data['x_Local_Migration']

        try:
            self.metadata[parser.sim_id] = parser.sim_data['x_Local_Migration']
        except:
            self.metadata[parser.sim_id] = 0
        # self.metadata[parser.sim_id] = parser.simulation.experiment.exp_name

        # parser.sim_data contains things like the collection_id's, sim_id's, environment (Belegost), etc.

        # print parser.raw_data[self.filenames[1]]["parameters"]["Simulation_Duration"]








    def finalize(self):
        # print self.my_data
        print ""

    def plot(self):
        import matplotlib.pyplot as plt
        import matplotlib.dates as mdates

        # start_date = "2001-01-01"  # Day 1 of simulation
        # date_format = "%Y-%m-%d"
        start_date = "2007-01-01"  # Day 1 of simulation
        date_format = "%Y-%m-%d"

        # Plot RDT prevalence as a function of time for node 0:
        # for sim_id, sim_data in self.my_data.items():
        #     plt.plot(sim_data)
        # for sim_id, data in self.prev_aggr.items():
        #     plt.plot(np.arange(self.n_tstep),self.prev_aggr[sim_id],label=self.metadata[sim_id]+' Prevalence')

        foo = mdates.strpdate2num(date_format)

        daynum = np.arange(self.n_tstep)
        daydates_list = []
        daydates_mdates = np.array([])
        for dayn in daynum:
            hold = convert_to_date(dayn, start_date, date_format=date_format)
            daydates_list.append(hold)
            daydates_mdates = np.append(daydates_mdates,foo(hold))
        # daydates2 = mdates.strpdate2num

        print daydates_mdates

        for sim_id, data in self.RDT_prev_aggr.items():
            # plt.plot(np.arange(self.n_tstep), self.RDT_prev_aggr[sim_id], label=self.metadata[sim_id]+' RDT Prevalence')

            # plt.plot(daynum, self.RDT_prev_aggr[sim_id], label=self.metadata[sim_id] + ' RDT Prevalence')
            plt.plot_date(daydates_mdates, self.RDT_prev_aggr[sim_id],fmt='-')#, label='x_Local_Migration = '+str(self.metadata[sim_id]),fmt='-')
        # plt.legend([s['environment'] for s in self.metadata.values()])




        # Plot Chiyabi prevalence data:
        if False:
            # Overplot Chiyabi RDT prevalence data:
            use_fulldate = True

            # Chiyabi Prevalence Data:
            # chiyabi_RDT_prev = {
            #     "1": {"date": "2012-06-15", "prev": 0.47242206},
            #     "2": {"date": "2012-08-21", "prev": 0.28082867},
            #     "3": {"date": "2012-10-31", "prev": 0.21063652},
            #     "4": {"date": "2013-06-08", "prev": 0.43310849},
            #     "5": {"date": "2013-08-14", "prev": 0.34187738},
            #     "6": {"date": "2013-07-28", "prev": 0.44410507},
            #     "7": {"date": "2014-12-11", "prev": 0.29185428},
            #     "8": {"date": "2015-01-31", "prev": 0.20071942},
            #     "9": {"date": "2015-09-09", "prev": 0.33360064},
            #     "10": {"date": "2015-12-11", "prev": 0.08664122}
            # }
            chiyabi_RDT_prev = {
                "1": {"date": "2012-07-01", "prev": 0.47242206},
                "2": {"date": "2012-09-30", "prev": 0.28082867},
                "3": {"date": "2012-11-30", "prev": 0.21063652},
                "4": {"date": "2013-07-01", "prev": 0.43310849},
                "5": {"date": "2013-08-31", "prev": 0.34187738},
                "6": {"date": "2013-10-31", "prev": 0.44410507},
                "7": {"date": "2014-12-31", "prev": 0.29185428},
                "8": {"date": "2015-03-01", "prev": 0.20071942},
                "9": {"date": "2015-09-30", "prev": 0.33360064},
                "10": {"date": "2016-02-29", "prev": 0.08664122}
            }


            chiyabi_RDT_prev_array = np.zeros([10,2])
            for round in chiyabi_RDT_prev.keys():
                # day_num = convert_to_day(chiyabi_RDT_prev[round]["date"], start_date, date_format=date_format)
                day_mdate = foo(chiyabi_RDT_prev[round]["date"])
                prev = chiyabi_RDT_prev[round]["prev"]
                # chiyabi_RDT_prev_array[int(round)-1] = np.array([day_num,prev])
                chiyabi_RDT_prev_array[int(round) - 1] = np.array([day_mdate, prev])

            plt.scatter(chiyabi_RDT_prev_array[:,0],chiyabi_RDT_prev_array[:,1],c='red',s=70,label='Data')






        # Plot Chiyabi interventions as vertical lines:
        if True:
            # plot_only_first = False

            # Event information files
            itn_event_file = "C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/cbever/chiyabi_hfca_itn_events.csv"
            irs_event_file = "C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/cbever/chiyabi_hfca_irs_events.csv"
            msat_event_file = "C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/cbever/chiyabi_hfca_msat_events.csv"
            mda_event_file = "C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/cbever/chiyabi_hfca_mda_events.csv"
            healthseek_event_file = "C:/Users/jsuresh//OneDrive - IDMOD/Code/zambia/cbever/chiyabi_hfca_healthseek_events.csv"
            stepd_event_file = "C:/Users/jsuresh//OneDrive - IDMOD/Code/zambia/cbever/chiyabi_hfca_stepd_events.csv"

            # Import event info
            itn_events = pd.read_csv(itn_event_file)
            irs_events = pd.read_csv(irs_event_file)
            msat_events = pd.read_csv(msat_event_file)
            mda_events = pd.read_csv(mda_event_file)
            healthseek_events = pd.read_csv(healthseek_event_file)
            stepd_events = pd.read_csv(stepd_event_file)

            # for date in itn_events['fulldate']: plt.axvline(convert_to_day(date, start_date, date_format=date_format), ls='dashed', color='C0') #,label='ITN Events')
            # for date in irs_events['fulldate']: plt.axvline(convert_to_day(date, start_date, date_format=date_format), ls='dashed', color='C1') #,label='IRS Events')
            # for date in msat_events['fulldate']: plt.axvline(convert_to_day(date, start_date, date_format=date_format), ls='dashed', color='C2') #,label='MSAT Events')
            # for date in mda_events['fulldate']: plt.axvline(convert_to_day(date, start_date, date_format=date_format), ls='dashed', color='C3') #,label='MDA Events')

            for date in itn_events['fulldate']: plt.axvline(foo(date), ls='dashed', color='C0') #,label='ITN Events')
            for date in irs_events['fulldate']: plt.axvline(foo(date), ls='dashed', color='C1') #,label='IRS Events')
            for date in msat_events['fulldate']:plt.axvline(foo(date), ls='dashed', color='C2') #,label='MSAT Events')
            for date in mda_events['fulldate']: plt.axvline(foo(date), ls='dashed', color='C3') #,label='MDA Events')


            # colors:
            # c0 = blue = ITN
            # c1 = green = IRS
            # c2 = red = MSAT
            # c3 = purple = MDA




        plt.legend()
        # plt.xlim([3000,7000])
        plt.xlim([foo("2010-01-01"), foo("2019-01-01")])
        plt.show()




        # Plot RDT prevalence for 9 largest nodes:
        if True:
            # Find 9 largest nodes:
            # pop_inits =
            pass


if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()

    am.add_experiment(retrieve_experiment("c4f143ce-d60c-e711-9400-f0921c16849c"))

    am.add_analyzer(RDTPrevAnalyzer())
    am.analyze()