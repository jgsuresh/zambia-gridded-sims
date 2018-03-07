from CalibMozambiqueAnalyzer import CalibMozambiqueAnalyzer
from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from relative_time import *
from simtools.AnalyzeManager.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
from CalibTestAnalyzer import CalibTestAnalyzer

if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()

    am.add_experiment(retrieve_experiment("3b7195e6-c020-e811-a2bf-c4346bcb7274"))


    am.add_analyzer(CalibMozambiqueAnalyzer("GriddedCalibSite"))
    am.analyze()