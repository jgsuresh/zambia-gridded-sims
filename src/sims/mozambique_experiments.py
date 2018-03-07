from experiment_setup import GriddedConfigBuilder

import numpy as np

from dtk.vector.species import set_species_param

from gridded_sim_general import *
from experiment_setup import CatchmentDemographicsGenerator

class MozambiqueExperiment(GriddedConfigBuilder):

    def __init__(self,
                 base,
                 exp_name,
                 catch,
                 healthseek_fn=None,
                 itn_fn=None,
                 irs_fn=None,
                 msat_fn=None,
                 mda_fn=None,
                 stepd_fn=None,
                 start_year=2001,
                 sim_length_years=19,
                 immunity_mode="naive",
                 num_cores=12,
                 parser_location='HPC'):

        self.catch = catch

        # Migration:
        # self.migration_on = True
        # self.gravity_migr_params = np.array([7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00])

        catch_cells = MozambiqueExperiment.find_cells_for_this_catchment(self.catch)

        super().__init__(base,
                         exp_name,
                         catch_cells,
                         region="Mozambique",
                         healthseek_fn=healthseek_fn,
                         itn_fn=itn_fn,
                         irs_fn=irs_fn,
                         msat_fn=msat_fn,
                         mda_fn=mda_fn,
                         stepd_fn=stepd_fn,
                         start_year=start_year,
                         sim_length_years=sim_length_years,
                         immunity_mode=immunity_mode,
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
            "TEMPORARY_RAINFALL": 2.2e7,
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
            "WATER_VEGETATION": 2e3,
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

    # def return_cb_for_calibration(self):
    #     self.implement_baseline_healthseeking()
    #     self.implement_interventions(self.cb, True, True, False, True, False)
    #     return self.cb

    def larval_params_func_for_calibration(self, grid_cells):
        return {"LINEAR_SPLINE": np.ones_like(grid_cells),
                "WATER_VEGETATION": np.ones_like(grid_cells),
                "TEMPORARY_RAINFALL": np.ones_like(grid_cells)}


    # Grid-cell/Node ID
    @staticmethod
    def find_cells_for_this_catchment(catch, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
        # Find which grid cells correspond to a given HFCA
        df = pd.read_csv(base + "data/mozambique/grid_lookup.csv")

        if catch == 'all':
            return np.array(df['grid_cell'])
        else:
            df_catch = df[df['catchment'] == catch]
            # df_catch = df[np.logical_or(df['catchment'] == catch.capitalize(),
            #                             df['catchment'] == catch.lower())]
            return np.array(df_catch['grid_cell'])

    @staticmethod
    def find_bairros_for_this_catchment(catch, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
        catch_cells = MozambiqueExperiment.find_cells_for_this_catchment(catch, base=base)

        df = pd.read_csv(base + "data/mozambique/grid_lookup_with_neighborhood.csv")
        in_catch = np.in1d(df['grid_cell'], catch_cells)

        bairro_name_list = sorted(list(set(df['bairro_name'][in_catch])))
        bairro_num_list = sorted(list(set(df['bairro'][in_catch])))
        num_bairros = len(bairro_num_list)

        dd = {"bairro_name_list": bairro_name_list,
              "bairro_num_list": bairro_num_list,
              "num_bairros": num_bairros}

        for bairro_num in bairro_num_list:
            in_bairro = df['bairro'] == bairro_num
            in_catch_bairro = np.logical_and(in_catch, in_bairro)
            dd[bairro_num] = sorted(list(set(df['grid_cell'][in_catch_bairro])))

        return dd


