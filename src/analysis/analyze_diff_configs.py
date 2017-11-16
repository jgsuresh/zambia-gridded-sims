from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from relative_time import *
from simtools.AnalyzeManager.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
import numpy as np
import pandas as pd

class ConfigDiffAnalyzer(BaseAnalyzer):

    filenames = ['config.json']

    def __init__(self):
        super(ConfigDiffAnalyzer, self).__init__()
        self.config_data = {}
        self.metadata = {}

    def filter(self, sim_metadata):
        return True

    def apply(self, parser):
        self.config_data[parser.sim_id] = parser.raw_data[self.filenames[0]]


    def finalize(self):
        pass

        sim_ids = self.config_data.keys()

        # d1 = self.config_data[sim_ids[0]]['parameters']
        # d2 = self.config_data[sim_ids[-1]]['parameters']
        d1 = self.config_data["e9533bd4-d60c-e711-9400-f0921c16849c"]['parameters']
        d2 = self.config_data["c1826603-4f87-e711-9401-f0921c16849d"]['parameters']


        self.dict_diff(d1,d2,sim_ids[0],sim_ids[1])

        # for key in d1:
        #     if key not in d2:
        #         print "Param {} in sim {} but not {}".format(key,sim_ids[0],sim_ids[1])
        #     else:
        #         if d1[key] == d2[key]:
        #             pass
        #         else:
        #
        #             print "Param {}: {} in {} ---- {} in {} ".format(key,d1[key],sim_ids[0],d2[key],sim_ids[1])
        #             print ""

        print ""

    def dict_diff(self,d1,d2,id1,id2):
        for key in d1:
            if key not in d2:
                print "Param {} in sim {} but not {}".format(key,id1,id2)
            else:
                if isinstance(d1[key],dict) and isinstance(d2[key],dict):
                    self.dict_diff(d1[key],d2[key],id1,id2)
                elif d1[key] == d2[key]:
                    pass
                else:
                    print "Param {}: {} in {} ---- {} in {} ".format(key,d1[key],id1,d2[key],id2)
                    print ""

    def plot(self):
        pass



if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()
    am.add_experiment(retrieve_experiment("c4f143ce-d60c-e711-9400-f0921c16849c"))
    am.add_experiment(retrieve_experiment("c0826603-4f87-e711-9401-f0921c16849d"))

    am.add_analyzer(ConfigDiffAnalyzer())
    am.analyze()