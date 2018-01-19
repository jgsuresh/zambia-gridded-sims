from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from relative_time import convert_to_day
from simtools.AnalyzeManager.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
import numpy as np
import pandas as pd

class MigrationAnalyzer(BaseAnalyzer):

    filenames = ['output\ReportEventRecorder.csv','output\SpatialReport_Population.bin']

    def __init__(self):
        super(MigrationAnalyzer, self).__init__()
        self.my_data = {}
        self.metadata = {}

        self.pop_init = {}
        self.n_travellers = {}
        self.n_trips = {}

    def filter(self, sim_metadata):
        # return sim_metadata['Run_Number'] == 0
        # print sim_metadata.keys()
        # print sim_metadata['run_number']
        return True
        # return parser.sim_data['Run_Number'] == 0

    def apply(self, parser):
        self.metadata[parser.sim_id] = parser.simulation.experiment.exp_name

        event_data = parser.raw_data[self.filenames[0]]
        pop_data = parser.raw_data[self.filenames[1]]

        # Get initial population of nodes:
        n_nodes = pop_data['n_nodes']

        pop_init_by_node = np.zeros(n_nodes)
        for ni in range(n_nodes):
            pop_init_by_node[ni] = pop_data['data'][0][ni]

        self.pop_init[parser.sim_id] = np.sum(pop_init_by_node)

        time_start = 1000
        time_window = 30

        event_data_windowed = event_data[np.logical_and(event_data['Time'] >= time_start, event_data['Time'] < time_start+time_window)]
        event_data_windowed = event_data_windowed[np.logical_or(event_data_windowed['Event_Name'] == 'Immigrating', event_data_windowed['Event_Name'] == 'Emigrating')]
        immig_events = event_data_windowed[event_data_windowed['Event_Name'] == 'Immigrating']

        traveller_ids = set(immig_events['Individual_ID'])
        self.n_travellers[parser.sim_id] = len(traveller_ids)

        self.n_trips[parser.sim_id] = np.zeros(self.n_travellers[parser.sim_id])

        i = 0
        for iid in traveller_ids:
            self.n_trips[parser.sim_id][i] = np.sum(immig_events['Individual_ID'] == iid)
            i += 1



    def finalize(self):
        # print self.my_data
        print("")

    def plot(self):
        import matplotlib.pyplot as plt

        # Plot histogram of trips
        for sim_id, data in self.n_trips.items():
            # data only contains data for travellers.  Need to add in "zero trips" for everyone who didn't travel.
            n_couch = self.pop_init[sim_id] - self.n_travellers[sim_id]
            full_data = np.append(data, np.zeros(int(n_couch)))
            plt.hist(full_data,histtype='stepfilled',alpha=0.4,log=True, label=self.metadata[sim_id])

        plt.legend()

        plt.show()



if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()

    am.add_experiment(retrieve_experiment("151f8b4b-867c-e711-9401-f0921c16849d"))

    am.add_analyzer(MigrationAnalyzer())
    am.analyze()