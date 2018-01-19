from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from relative_time import *
from simtools.AnalyzeManager.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
import numpy as np
import pandas as pd

class AgeStratificationAnalyzer(BaseAnalyzer):

    filenames = ['output\MalariaSummaryReport_AnnualAverage.json']



    def __init__(self):
        super(AgeStratificationAnalyzer, self).__init__()
        self.report_times = None
        self.age_bins = None
        self.pop_data = {}
        self.raw_pop_data = {}
        self.tot_pop = {}
        self.metadata = {}

        # self.plot_labels = ["Multinode x_Local_Migration=10","Multinode x_Local_Migration=1","Multinode x_Local_Migration=0.1","Multinode x_Local_Migration=0", "Singlenode"]

    def filter(self, sim_metadata):
        return sim_metadata['Run_Number'] == 0
        # if isinstance(sim_metadata['IRS'],bool):
        #     all_intervene_on = (sim_metadata['IRS'] and sim_metadata['ITNs'] and sim_metadata['MDA'] and sim_metadata['MSAT'] and sim_metadata['StepD'] and sim_metadata['Healthseek'])
        # else:
        #     all_intervene_on = (sim_metadata['IRS'] == "true" and
        #                      sim_metadata['ITNs']  == "true" and
        #                      sim_metadata['MDA']  == "true" and
        #                      sim_metadata['MSAT']  == "true" and
        #                      sim_metadata['StepD']  == "true" and
        #                      sim_metadata['Healthseek'] == "true")
        #
        #
        # return all_intervene_on


    def apply(self, parser):
        raw = parser.raw_data[self.filenames[0]]
        exp_id = parser.experiment.exp_id

        self.report_times = raw['DataByTime']['Time Of Report']
        self.age_bins = raw['Metadata']['Age Bins']
        raw_pop_data= raw['DataByTimeAndAgeBins']['Average Population by Age Bin']
        self.raw_pop_data[exp_id] = raw_pop_data

        self.pop_data[exp_id] = {}
        for ri in range(len(self.report_times)):
            for ai in range(len(self.age_bins)):
                if ri == 0:
                    self.pop_data[exp_id][ai] = np.array([raw_pop_data[ri][ai]])
                else:
                    self.pop_data[exp_id][ai] = np.append(self.pop_data[exp_id][ai],raw_pop_data[ri][ai])

        # Sanity check:
        self.tot_pop[exp_id] = np.zeros(0)
        for ri in range(len(self.report_times)):
            pop_this_timestep = np.array(raw_pop_data[ri])
            print("Total population at report time {} = {}".format(ri,np.sum(pop_this_timestep)))
            self.tot_pop[exp_id] = np.append(self.tot_pop[exp_id],np.sum(pop_this_timestep))



    def finalize(self):
        print("")

    def plot(self):
        import matplotlib.pyplot as plt
        import seaborn
        import matplotlib.dates as mdates


        if False:
            plt.close('all')

            ax=plt.subplot(2,2,1)
            ax.set_title("Age 0-4")
            for exp_id in self.pop_data.keys():
                if exp_id == self.pop_data.keys()[0]:
                    ls = '-'
                    label_flag = 1
                else:
                    ls = '--'
                    label_flag = 0

                for ai in range(len(self.age_bins)):
                    if label_flag == 1:
                        label = self.age_bins[ai]
                    else:
                        label = None

                    if ai >= 0 and ai < 6:
                        plt.plot(self.report_times[:-1], self.pop_data[exp_id][ai][:-1]/self.tot_pop[exp_id][:-1],linestyle=ls, label=label, alpha=0.7)
            plt.legend()

            ax =plt.subplot(2,2,2)
            ax.set_title("Age 4-10")
            for exp_id in self.pop_data.keys():
                if exp_id == self.pop_data.keys()[0]:
                    ls = '-'
                    label_flag = 1
                else:
                    ls = '--'
                    label_flag = 0

                for ai in range(len(self.age_bins)):
                    if label_flag == 1:
                        label = self.age_bins[ai]
                    else:
                        label = None

                    if ai >= 6 and ai < 12:
                        plt.plot(self.report_times[:-1], self.pop_data[exp_id][ai][:-1]/self.tot_pop[exp_id][:-1],linestyle=ls, label=label, alpha=0.7)
            plt.legend()

            ax =plt.subplot(2,2,3)
            ax.set_title("Age 10-20")
            for exp_id in self.pop_data.keys():
                if exp_id == self.pop_data.keys()[0]:
                    ls = '-'
                    label_flag = 1
                else:
                    ls = '--'
                    label_flag = 0

                for ai in range(len(self.age_bins)):
                    if label_flag == 1:
                        label = self.age_bins[ai]
                    else:
                        label = None

                    if ai >= 12 and ai < 18:
                        plt.plot(self.report_times[:-1], self.pop_data[exp_id][ai][:-1]/self.tot_pop[exp_id][:-1],linestyle=ls, label=label, alpha=0.7)
            plt.legend()

            ax =plt.subplot(2,2,4)
            ax.set_title("Age 20-1000")
            for exp_id in self.pop_data.keys():
                if exp_id == self.pop_data.keys()[0]:
                    ls = '-'
                    label_flag = 1
                else:
                    ls = '--'
                    label_flag = 0

                for ai in range(len(self.age_bins)):
                    if label_flag == 1:
                        label = self.age_bins[ai]
                    else:
                        label = None

                    if ai >= 18 and ai < 24:
                        plt.plot(self.report_times[:-1], self.pop_data[exp_id][ai][:-1]/self.tot_pop[exp_id][:-1],linestyle=ls, label=label, alpha=0.7)
            plt.legend()

            plt.show()


        if True:
            # Plot age distribution at early and late timepoints
            ax=plt.subplot(2,1,1)
            for exp_id in self.pop_data.keys():
                plt.plot(np.array(self.raw_pop_data[exp_id][0])/self.tot_pop[exp_id][0])
            ax.set_title("Early")
            ax.set_xticks(range(24))
            ax.set_xticklabels(self.age_bins)

            ax=plt.subplot(2,1,2)
            for exp_id in self.pop_data.keys():
                plt.plot(np.array(self.raw_pop_data[exp_id][-2])/self.tot_pop[exp_id][-2])
            ax.set_title("Late")
            ax.set_xticks(range(24))
            ax.set_xticklabels(self.age_bins)
            plt.show()


        #
        # for exp_id in self.pop_data.keys():
        #     plt.plot_date(self.report_times, self.pop_data[exp_id],fmt='-',c=c,linewidth=lw,label=label,alpha=0.4)
        #     plt.plot_date(self.report_times, self.pop_data[exp_id], fmt='-', c=c, linewidth=lw, label=label, alpha=0.4)
        #     plt.plot_date(self.report_times, self.pop_data[exp_id], fmt='-', c=c, linewidth=lw, label=label, alpha=0.4)
        #     plt.plot_date(self.report_times, self.pop_data[exp_id], fmt='-', c=c, linewidth=lw, label=label, alpha=0.4)
        # plt.legend([s['environment'] for s in self.metadata.values()])






if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()

    am.add_experiment(retrieve_experiment("f4ecdcc6-768c-e711-9401-f0921c16849d")) # L1
    # am.add_experiment(retrieve_experiment("001a9f44-758c-e711-9401-f0921c16849d")) # L5
    am.add_experiment(retrieve_experiment("4188b9de-e28c-e711-9401-f0921c16849d"))  # L6

    am.add_analyzer(AgeStratificationAnalyzer())
    am.analyze()