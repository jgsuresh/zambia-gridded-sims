from calibtool.CalibSite import CalibSite
from GriddedRDTLikelihoodAnalyzer import GriddedRDTLikelihoodAnalyzer
from CalibTestAnalyzer import CalibTestAnalyzer

import os
import pandas as pd

class GriddedCalibSite(CalibSite):
    def __init__(self):
        super().__init__("GriddedCalibSite")

    # def get_reference_data(self, reference_type):
    #     # Load the RDT data
    #     #todo Make this more general, for other users
    #     base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
    #     reference_csv = base + "data/interventions/kariba/2018-01-23/raw/grid_prevalence_with_dates.csv"
    #     reference_data = pd.read_csv(reference_csv)
    #     return reference_data

    def get_analyzers(self):
        # return GriddedRDTLikelihoodAnalyzer(self),
        return CalibTestAnalyzer(self),

    # def _load_RDT_prev_data(self, reference_csv):
    #     reference_data_df = pd.read_csv(reference_csv)
    #
    #     reference_data_df
