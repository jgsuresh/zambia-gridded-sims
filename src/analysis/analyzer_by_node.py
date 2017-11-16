from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from relative_time import *
from simtools.AnalyzeManager.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
import numpy as np
import pandas as pd

class NodeAnalyzer(BaseAnalyzer):

    filenames = ['output\SpatialReport_Population.bin', 'output\SpatialReport_Prevalence.bin', 'output\SpatialReport_New_Diagnostic_Prevalence.bin', 'output\SpatialReport_Daily_Bites_Per_Human.bin']
    # filenames = ['output\InsetChart.json']

    def __init__(self):
        super(NodeAnalyzer, self).__init__()
        self.bite_data = {}
        self.RDT_prev_data = {}
        self.metadata = {}

    def filter(self, sim_metadata):
        if isinstance(sim_metadata['IRS'],bool):
            all_intervene_on = (sim_metadata['IRS'] and sim_metadata['ITNs'] and sim_metadata['MDA'] and sim_metadata['MSAT'] and sim_metadata['StepD'] and sim_metadata['Healthseek'])
        else:
            all_intervene_on = (sim_metadata['IRS'] == "true" and
                             sim_metadata['ITNs']  == "true" and
                             sim_metadata['MDA']  == "true" and
                             sim_metadata['MSAT']  == "true" and
                             sim_metadata['StepD']  == "true" and
                             sim_metadata['Healthseek'] == "true")


        return all_intervene_on


    def apply(self, parser):
        pop_raw = parser.raw_data[self.filenames[0]]
        prev_raw = parser.raw_data[self.filenames[1]]
        RDT_prev_raw = parser.raw_data[self.filenames[2]]
        bite_raw = parser.raw_data[self.filenames[3]]
        exp_id = parser.experiment.exp_id
        sim_id = parser.sim_id

        self.n_tstep = pop_raw['n_tstep']
        self.n_nodes = pop_raw['n_nodes']
        node_ids = pop_raw['nodeids']

        # Get initial population of nodes:
        self.pop_init = np.zeros(self.n_nodes)
        for ni in range(self.n_nodes):
            self.pop_init[ni] = pop_raw['data'][0][ni]

        # Ensure that each exp_id is in the data dictionaries
        for data in [self.RDT_prev_data,self.bite_data]:
            if exp_id not in data:
                data[exp_id] = {}



            # self.prev_aggr[parser.sim_id][i] = np.sum(pop_data['data'][i]*prev_data['data'][i])/np.sum(pop_data['data'][i])
            # self.RDT_prev_aggr[parser.sim_id][i] = np.sum(pop_data['data'][i] * RDT_prev_data['data'][i]) / np.sum(pop_data['data'][i])

        pass
        # foo = pd.DataFrame.from_dict(bite_raw['data'],orient='index')
        # print foo.columns

        print "done"
        # for node_id in node_ids:
        #     bite_arr = np.zeros(self.n_tstep)
        #     RDT_prev_arr = np.zeros(self.n_tstep)
        #     for ti in np.arange(self.n_tstep):
        #         bite_arr[ti] = 0
        # self.bite_data[exp_id][sim_id] = bite_raw['Channels']['Daily Bites per Human']['Data']
        # self.RDT_prev_data[exp_id][sim_id] = raw['Channels']['New Diagnostic Prevalence']['Data']
        #
        # # self.bite_data[parser.sim_id] = raw['Channels']['Daily Bites per Human']['Data']
        # # self.RDT_prev_data[parser.sim_id] = raw['Channels']['New Diagnostic Prevalence']['Data']
        #
        # self.n_tstep = np.size(self.bite_data[exp_id][parser.sim_id])


    def finalize(self):
        print ""

    def plot(self):
        import matplotlib.pyplot as plt
        import seaborn
        import matplotlib.dates as mdates

        plt.close('all')

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
        # daydates2 = mdates.strpdate2num

        print daydates_mdates

        jflag = 0
        jflag2 = 0

        for exp_id in self.RDT_prev_data.keys():
            for sim_id, data in self.RDT_prev_data[exp_id].items():
                if sim_id == "e9533bd4-d60c-e711-9400-f0921c16849c":
                    c = 'C5'
                    lw =1.2
                    label='Caitlin'
                elif exp_id == "16a8495d-1d89-e711-9401-f0921c16849d":
                    c = 'C1'
                    lw = 0.2
                    label=None
                    if jflag2 == 0:
                        jflag2 = 1
                        label='old .exe'
                else:
                    c = 'black'
                    lw = 1.2
                    label = None
                    if jflag == 0:
                        jflag = 1
                        label = 'new .exe' #'Josh'
                plt.plot_date(daydates_mdates, self.RDT_prev_data[exp_id][sim_id],fmt='-',c=c,linewidth=lw,label=label)
        plt.legend([s['environment'] for s in self.metadata.values()])




        # Plot Chiyabi prevalence data:
        if True:
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

            label_flag = 0
            for date in itn_events['fulldate']:
                label = None
                if label_flag == 0:
                    label = "ITNs"
                    label_flag = 1
                plt.axvline(foo(date), ls='dashed', color='C0',label=label)

            label_flag = 0
            for date in irs_events['fulldate']:
                label = None
                if label_flag == 0:
                    label = "IRS"
                    label_flag = 1
                plt.axvline(foo(date), ls='dashed', color='C1', label=label)

            label_flag = 0
            for date in msat_events['fulldate']:
                label = None
                if label_flag == 0:
                    label = "MSAT"
                    label_flag = 1
                plt.axvline(foo(date), ls='dashed', color='C2', label=label)

            label_flag = 0
            for date in mda_events['fulldate']:
                label = None
                if label_flag == 0:
                    label = "MDA"
                    label_flag = 1
                plt.axvline(foo(date), ls='dashed', color='C3', label=label)
            # for date in irs_events['fulldate']: plt.axvline(foo(date), ls='dashed', color='C1') #,label='IRS Events')
            # for date in msat_events['fulldate']:plt.axvline(foo(date), ls='dashed', color='C2') #,label='MSAT Events')
            # for date in mda_events['fulldate']: plt.axvline(foo(date), ls='dashed', color='C3') #,label='MDA Events')


            # colors:
            # c0 = blue = ITN
            # c1 = green = IRS
            # c2 = red = MSAT
            # c3 = purple = MDA




        plt.legend()
        # plt.xlim([3000,7000])
        plt.title("Chiyabi (Single-Node)")
        plt.xlabel("Date")
        plt.ylabel("RDT Prevalence")
        plt.xlim([foo("2010-01-01"), foo("2019-01-01")])
        plt.show()




if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()

    am.add_experiment(retrieve_experiment("96afb603-1789-e711-9401-f0921c16849d")) #Josh multinode
    am.add_analyzer(NodeAnalyzer())
    am.analyze()