from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from relative_time import *
from simtools.AnalyzeManager.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
import numpy as np
import pandas as pd

class RDTPrevAnalyzer(BaseAnalyzer):

    # filenames = ['output\SpatialReport_Population.bin', 'output\SpatialReport_Prevalence.bin', 'output\SpatialReport_New_Diagnostic_Prevalence.bin']
    filenames = ['output\InsetChart.json']

    def __init__(self):
        super(RDTPrevAnalyzer, self).__init__()
        self.bite_data = {}
        self.RDT_prev_data = {}
        self.metadata = {}

    def filter(self, sim_metadata):
        # return sim_metadata['Run_Number'] == 0
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
        raw = parser.raw_data[self.filenames[0]]
        exp_id = parser.experiment.exp_id

        if exp_id not in self.bite_data:
            self.bite_data[exp_id] = {}
            self.RDT_prev_data[exp_id] = {}

        self.bite_data[exp_id][parser.sim_id] = raw['Channels']['Daily Bites per Human']['Data']
        self.RDT_prev_data[exp_id][parser.sim_id] = raw['Channels']['New Diagnostic Prevalence']['Data']

        self.n_tstep = np.size(self.bite_data[exp_id][parser.sim_id])


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


        for sp in [0,1]:
            if sp == 0:
                data_dict = self.RDT_prev_data
                plot_chiyabi_rdt = True
                ax = plt.subplot(2,1,1)
            elif sp == 1:
                data_dict = self.bite_data
                plot_chiyabi_rdt = False
                ax = plt.subplot(2,1,2)


            jflag = 0
            jflag2 = 0


            for exp_id in data_dict.keys():
                for sim_id, data in data_dict[exp_id].items():
                    if sim_id == "e9533bd4-d60c-e711-9400-f0921c16849c":
                        c = 'C5'
                        lw =1.5
                        label='Caitlin'
                        alpha=1.0
                    else:
                        c = 'black'
                        lw = 1.2
                        label = None
                        alpha=0.4
                        if jflag == 0:
                            jflag = 1
                            label = 'Josh'
                    plt.plot_date(daydates_mdates, data_dict[exp_id][sim_id],fmt='-',c=c,linewidth=lw,label=label,alpha=alpha)


            # Plot Chiyabi prevalence data:
            if plot_chiyabi_rdt:

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

                # colors:
                # c0 = blue = ITN
                # c1 = green = IRS
                # c2 = red = MSAT
                # c3 = purple = MDA




            ax.set_xlabel("Date")
            if sp == 0:
                ax.set_ylabel("RDT Prevalence")
                ax.set_title("Chiyabi (Single-Node)")
            elif sp == 1:
                ax.set_ylabel("Daily Bites per Human")
                plt.legend(loc=0)

            ax.set_xlim([foo("2010-01-01"), foo("2019-01-01")])

        plt.show()




if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()

    am.add_experiment(retrieve_experiment("c4f143ce-d60c-e711-9400-f0921c16849c")) #Cait
    # am.add_experiment(retrieve_experiment("40221d78-b289-e711-9401-f0921c16849d"))  # old .exe with new climate
    # am.add_experiment(retrieve_experiment("5e854ab1-bb89-e711-9401-f0921c16849d")) # old .exe with old climate

    # am.add_experiment(retrieve_experiment("40221d78-b289-e711-9401-f0921c16849d"))  # old .exe with old climate and enable_demographics_other = 0
    # am.add_experiment(retrieve_experiment("5e854ab1-bb89-e711-9401-f0921c16849d")) # old .exe with old climate and enable_demographics_other = 1


    # am.add_experiment(retrieve_experiment("bc6d02eb-e789-e711-9401-f0921c16849d")) #cait_rerun
    # am.add_experiment(retrieve_experiment("14584d58-ee89-e711-9401-f0921c16849d"))  # cait_rerun_newdemo
    # am.add_experiment(retrieve_experiment("c27fe675-f089-e711-9401-f0921c16849d"))  # cait_rerun_newdemo_unfancyweather

    # am.add_experiment(retrieve_experiment("4a2a24b5-f589-e711-9401-f0921c16849d")) # deconstruct_L1
    # am.add_experiment(retrieve_experiment("0e95c4e4-f589-e711-9401-f0921c16849d"))  # deconstruct_L2
    # am.add_experiment(retrieve_experiment("316f7813-f689-e711-9401-f0921c16849d"))  # deconstruct_L3
    # am.add_experiment(retrieve_experiment("7346abff-f689-e711-9401-f0921c16849d"))  # deconstruct_L4
    # am.add_experiment(retrieve_experiment("f2e1a7e7-308c-e711-9401-f0921c16849d"))  # deconstruct_L4

    # am.add_experiment(retrieve_experiment("7346abff-f689-e711-9401-f0921c16849d")) # L4
    # am.add_experiment(retrieve_experiment("001a9f44-758c-e711-9401-f0921c16849d")) #L5

    # am.add_experiment(retrieve_experiment("2007274b-3c8d-e711-9401-f0921c16849d")) # L10 4761
    # am.add_experiment(retrieve_experiment("fd452d45-3c8d-e711-9401-f0921c16849d"))  # L10 1000

    am.add_experiment(retrieve_experiment("ed7c56e7-f58d-e711-9401-f0921c16849d"))  # L10 4761 NEW INFECTIVITY

    am.add_analyzer(RDTPrevAnalyzer())
    am.analyze()