from experiment_setup import COMPS_Experiment

import numpy as np

from dtk.vector.species import set_species_param

class Zambia_Experiment(COMPS_Experiment):
    def __init__(self,
                 base,
                 exp_name,
                 grid_pop_csv_file,
                 start_year=2001,
                 sim_length_years=19,
                 num_cores=12,
                 parser_location='HPC'):

        # Migration:
        self.migration_on = True
        self.gravity_migr_params = np.array([7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00])

        super().__init__(base,
                         exp_name,
                         grid_pop_csv_file,
                         start_year=start_year,
                         sim_length_years=sim_length_years,
                         num_cores=num_cores,
                         parser_location=parser_location)

        self.mozambique_setup()

        # Migration amplitude:
        self.cb.update_params({
            "x_Local_Migration": 4
        })

    def mozambique_setup(self):
        # Uses vector splines from Prashanth's Mozambique entomology calibration

        self.africa_setup()

        # Vector properties:
        self.cb.update_params({'Vector_Species_Names': ['arabiensis', 'funestus']})

        # Arabiensis
        set_species_param(self.cb, 'arabiensis', 'Larval_Habitat_Types', {
            "LINEAR_SPLINE": {
                "Capacity_Distribution_Per_Year": {
                    "Times": [0.0, 30.417, 60.833, 91.25, 121.667, 152.083, 182.5, 212.917, 243.333, 273.75, 304.167,
                              334.583],
                    "Values": [0.273953355,
                               4.226011848,
                               5.140191814,
                               9.363408701,
                               0.0,
                               0.414082115,
                               0.139915067,
                               0.186456901,
                               0.015611024,
                               0.101027567,
                               0.0,
                               0.121014426
                               ]
                },
                "Max_Larval_Capacity": pow(10, 8.5) # Uncertain-- this exponent was what we calibrate for
            }
        })

        # Funestus
        set_species_param(self.cb, 'funestus', 'Larval_Habitat_Types', {
            "LINEAR_SPLINE": {
                "Capacity_Distribution_Per_Year": {
                    "Times": [
                        0.0,
                        30.417,
                        60.833,
                        91.25,
                        121.667,
                        152.083,
                        182.5,
                        212.917,
                        243.333,
                        273.75,
                        304.167,
                        334.583
                    ],
                    "Values": [0.0,
                               1.202730029,
                               0.112447779,
                               1.467850365,
                               2.470962168,
                               1.064668156,
                               4.806719314,
                               0.914212162,
                               9.919572963,
                               0.437353893,
                               0.392657387,
                               1.213697659
                               ]
                },
                "Max_Larval_Capacity": pow(10, 8.5) # Uncertain-- this exponent was what we calibrate for
            }
        })

    def return_cb_for_calibration(self):
        self.implement_baseline_healthseeking()
        self.implement_interventions(self.cb, True, True, False, True, False)
        return self.cb