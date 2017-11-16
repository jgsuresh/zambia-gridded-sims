from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from simtools.AnalyzeManager.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment


class EIRAnalyzer(BaseAnalyzer):

    filenames = ['output\\SpatialReport_Daily_EIR.bin', 'config.json']

    def __init__(self):
        super(EIRAnalyzer, self).__init__()
        self.my_data = {}
        self.metadata = {}

    def filter(self, sim_metadata):
        return sim_metadata['environment'] == "Belegost"

    def apply(self, parser):
        data = parser.raw_data[self.filenames[0]]
        self.my_data[parser.sim_id] = data['data'][500]
        self.metadata[parser.sim_id] = parser.sim_data

        print parser.raw_data[self.filenames[1]]["parameters"]["Simulation_Duration"]

    def finalize(self):
        # print self.my_data
        print ""

    def plot(self):
        import matplotlib.pyplot as plt
        for sim_id, sim_data in self.my_data.items():
            plt.plot(sim_data)

        plt.legend([s['environment'] for s in self.metadata.values()])
        plt.show()




if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()
    am.add_experiment(retrieve_experiment('0b9a7d63-8471-e711-9401-f0921c16849d'))
    am.add_experiment(retrieve_experiment('a08bf1da-b670-e711-9401-f0921c16849d'))
    am.add_analyzer(EIRAnalyzer())
    am.analyze()