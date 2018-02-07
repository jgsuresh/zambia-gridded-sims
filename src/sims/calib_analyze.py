from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from relative_time import *
from simtools.AnalyzeManager.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
from CalibTestAnalyzer import CalibTestAnalyzer

if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()

    am.add_experiment(retrieve_experiment("0c8223c6-7a0b-e811-9415-f0921c16b9e5"))


    am.add_analyzer(CalibTestAnalyzer("GriddedCalibSite"))
    am.analyze()