from experiment_setup import GriddedConfigBuilder

import numpy as np

from dtk.vector.species import set_species_param

from gridded_sim_general import *
from experiment_setup import CatchmentDemographicsGenerator

class ZambiaExperiment(GriddedConfigBuilder):

    base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'

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
        self.region = "Zambia"

        # Migration:
        # self.migration_on = True
        # self.gravity_migr_params = np.array([7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00])

        catch_cells = ZambiaExperiment.find_cells_for_this_catchment(self.catch)

        super().__init__(base,
                         exp_name,
                         catch_cells,
                         region="Zambia",
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

        self.zambia_setup()

        # Migration amplitude:
        self.cb.update_params({
            "x_Local_Migration": 4
        })


    def zambia_setup(self):
        self.africa_setup()

        # Vector properties:
        self.cb.update_params({'Vector_Species_Names': ['arabiensis', 'funestus']})

        # Arabiensis
        set_species_param(self.cb, 'arabiensis', 'Larval_Habitat_Types', {
            "CONSTANT": 2000000.0,
            "TEMPORARY_RAINFALL": 100000000.0
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
                    "Values": [
                        0.0,
                        0.0,
                        0.0,
                        0.2,
                        0.8,
                        1.0,
                        1.0,
                        1.0,
                        0.5,
                        0.2,
                        0.0,
                        0.0
                    ]
                },
                "Max_Larval_Capacity": 100000000.0
            },
            "WATER_VEGETATION": 2000000.0
        })

    def gen_immunity_file_from_milen_clusters(self): #FIXME
        # Inputs:
        #   --demographics file with Caitlin-grid cell IDs as the node labels
        #   --CSV which maps from Caitlin-grid cell IDs to Milen-cluster IDs
        #   --CSV which maps from Milen-cluster IDs to climate category (e.g. Sinamalima, Gwembe, etc.)
        #   --Larval habitat parameters for same Milen-cluster IDs
        #   --Place to find Milen's immunity files (so we can use the one that's the closest match for the given cluster)

        # Outputs:
        #   --multi-node immunity file where each node's immunity is set to the best-fit approximation by Milen's calibration

        def generate_immun_dict_from_demo_file(cell_ids, node_ids):
            n_nodes = len(cell_ids)
            immun_fn_list = ZambiaExperiment.closest_milen_immunity_overlay_filenames_for_grid_cells(cell_ids)

            d = {}
            d["Nodes"] = []
            d["Defaults"] = {}
            d["Metadata"] = {}
            d["Metadata"]["Author"] = "Josh Suresh"
            d["Metadata"]["IdReference"] = "Gridded world grump30arcsec"
            d["Metadata"]["NodeCount"] = n_nodes

            for i in range(n_nodes):
                immun_fn = immun_fn_list[i]
                f = open(immun_fn, 'r')
                imm_dict = json.load(f)
                f.close()

                # node = {}
                # node["NodeAttributes"] = imm_dict['Defaults'].copy()
                node = imm_dict['Defaults'].copy()
                node["NodeID"] = node_ids[i]
                d["Nodes"].append(node)

            return d


        # Open demographics file
        with open(self.demo_fp_full, 'r') as f:
            demo_dict = json.load(f)

        # Get NodeID and grid cell (encoded as FacilityName)
        node_ids = []
        grid_cell_ids = []
        for node in demo_dict['Nodes']:
            node_ids.append(node['NodeID'])
            grid_cell_ids.append(int(node['NodeAttributes']['FacilityName']))

        # Get filenames of immunity file which correspond to best-fit of larval params and climate category for this grid cell:
        imm_dict = generate_immun_dict_from_demo_file(grid_cell_ids, node_ids)

        # Dump as new json file
        with open(self.immun_fp_full, 'w') as f:
            json.dump(imm_dict, f, indent=4)

    def get_immunity_file_from_single_node(self):
        # Inputs:
        #   --demographics file
        #   --single-node immunity file (e.g. "Sinamalima_1_node_immune_init_p1_33_p2_117.json")

        # Outputs:
        #   --multi-node immunity file where each node has identical immunity values.

        # Open demographics file
        with open(self.demo_fp_full, 'r') as f:
            demo_dict = json.load(f)

        # Get NodeID list
        node_ids = []
        for node in demo_dict['Nodes']:
            node_ids.append(node['NodeID'])

        # Open single-node immunity file
        with open(self.imm_1node_fp, 'r') as f:
            imm_dict = json.load(f)

        # Edit node list in this dictionary
        imm_node_list = imm_dict['Nodes']
        for node_id in node_ids:
            imm_node_list.append({u'NodeID': node_id})

        del imm_node_list[0]  # Remove the original node that was in the single-node immunity file

        # Edit node metadata to reflect new number of nodes:
        imm_dict['Metadata']['NodeCount'] = len(imm_node_list)

        # Dump as new json file
        with open(self.immun_fp_full, 'w') as f:
            json.dump(imm_dict, f, indent=4)

    # def add_larval_habitats_to_demo(self, demo_dict, larval_params=None):
    #     # Add larval habitat multipliers to demographics file
    #
    #     def add_larval_habitat_to_node(node_item, const_h, temp_h, water_h, linear_h):
    #         calib_single_node_pop = 1000
    #
    #         # This is now done in the demographics generator itself:
    #         # birth_rate = (float(node_item['NodeAttributes']['InitialPopulation']) / (1000 + 0.0)) * 0.12329
    #         # node_item['NodeAttributes']['BirthRate'] = birth_rate
    #
    #         pop_multiplier = float(node_item['NodeAttributes']['InitialPopulation']) / (calib_single_node_pop + 0.0)
    #
    #         temp_multiplier = temp_h * pop_multiplier
    #         linear_multiplier = linear_h * pop_multiplier
    #         const_multiplier = const_h * pop_multiplier  # NOTE: Pre 1/24/2018, this was not set to scale with population
    #         water_multiplier = water_h * pop_multiplier
    #
    #         node_item['NodeAttributes']['LarvalHabitatMultiplier'] = {
    #             "CONSTANT": const_multiplier,
    #             "TEMPORARY_RAINFALL": temp_multiplier,
    #             "WATER_VEGETATION": water_multiplier,
    #             "LINEAR_SPLINE": linear_multiplier
    #         }
    #
    #     if self.larval_params_mode == 'milen':
    #         # get grid cells from pop csv file:
    #         # for those grid cells, get corresponding arab/funest params
    #         # loop over nodes [order will correspond, by construction, to pop csv ordering]
    #         # give each node the corresponding larval params
    #
    #         # Load pop csv file to get grid cell numbers:
    #         # pop_df = pd.read_csv(self.grid_pop_csv_file)
    #         # grid_cells = np.array(pop_df['node_label'])
    #         catch_cells = ZambiaExperiment.find_cells_for_this_catchment(self.catch, path_from_base="data/mozambique/grid_lookup.csv")
    #
    #         # From those grid cells, and the Milen-clusters they correspond to, get best-fit larval habitat parameters
    #         arab_params, funest_params = ZambiaExperiment.find_milen_larval_param_fit_for_grid_cells(catch_cells)
    #
    #         # Loop over nodes in demographics file (which will, by construction, correspond to the grid pop csv ordering)
    #         i = 0
    #         for node_item in demo_dict['Nodes']:
    #             # if larval_params:
    #             #     const_h = larval_params['const_h']
    #             #     temp_h = larval_params['temp_h']
    #             #     water_h = larval_params['water_h']
    #             #     linear_h = larval_params['linear_h']
    #             # else:
    #             const_h = 1.
    #             temp_h = arab_params[i]
    #             water_h = 1.
    #             linear_h = funest_params[i]
    #
    #             add_larval_habitat_to_node(node_item, const_h, temp_h, water_h, linear_h)
    #
    #             i += 1
    #
    #     elif self.larval_params_mode == 'uniform':
    #         if larval_params:
    #             const_h = larval_params['const_h']
    #             temp_h = larval_params['temp_h']
    #             water_h = larval_params['water_h']
    #             linear_h = larval_params['linear_h']
    #         else:
    #             const_h = 1.
    #             temp_h = 122.
    #             water_h = 1.
    #             linear_h = 97.
    #
    #         for node_item in demo_dict['Nodes']:
    #             add_larval_habitat_to_node(node_item, const_h, temp_h, water_h, linear_h)
    #
    #
    #     elif self.larval_params_mode == 'calibrate':
    #         for node_item in demo_dict['Nodes']:
    #             add_larval_habitat_to_node(node_item, 1, 1, 1, 1)
    #
    #     return demo_dict

    def larval_params_func_milen(self, grid_cells):
        # get grid cells from pop csv file:
        # for those grid cells, get corresponding arab/funest params
        # loop over nodes [order will correspond, by construction, to pop csv ordering]
        # give each node the corresponding larval params

        # Load pop csv file to get grid cell numbers:
        # pop_df = pd.read_csv(self.grid_pop_csv_file)
        # grid_cells = np.array(pop_df['node_label'])
        # From those grid cells, and the Milen-clusters they correspond to, get best-fit larval habitat parameters
        arab_params, funest_params = ZambiaExperiment.find_milen_larval_param_fit_for_grid_cells(grid_cells)

        return {"CONSTANT": np.ones_like(grid_cells),
                "TEMPORARY_RAINFALL": arab_params,
                "LINEAR_SPLINE": funest_params,
                "WATER_VEGETATION": np.ones_like(grid_cells)}

    # Grid-cell/Node ID
    @staticmethod
    def find_cells_for_this_catchment(catch, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
        # Find which grid cells correspond to a given HFCA
        df = pd.read_csv(base + "data/zambia/grid_lookup.csv")

        if catch == 'all':
            return np.array(df['grid_cell'])
        else:
            df_catch = df[df['catchment'] == catch]
            # df_catch = df[np.logical_or(df['catchment'] == catch.capitalize(),
            #                             df['catchment'] == catch.lower())]
            return np.array(df_catch['grid_cell'])

    # Milen clusters
    @staticmethod
    def find_milen_cluster_for_grid_cells(cell_ids, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
        # Load Caitlin's cluster_to_grid_lookup.csv to identify, for each grid cell, which Milen-cluster it belongs to
        # NOTE: Assumes that cells is sorted
        cl_lookup = pd.read_csv(base + "data/zambia/milen_clusters/cluster_to_grid_lookup.csv")

        # fixme This assumes that every possible cell in cell_ids is in the lookup_cells.  Fails if not the case.
        # lookup_cells = np.array(cl_lookup['grid_cell'])
        # lookup_id = np.array(cl_lookup['cluster_id'])
        # in_desired_cells = np.in1d(lookup_cells,cell_ids)
        # cluster_ids = list(lookup_id[in_desired_cells])

        hold_df = pd.DataFrame({
            "cell_ids": cell_ids
        })

        hold_df = hold_df.merge(cl_lookup, how='left', left_on='cell_ids', right_on='grid_cell')
        cluster_ids = hold_df['cluster_id']
        cluster_ids = cluster_ids.fillna(method='ffill')  # Interpolate for any missing cells
        cluster_ids = cluster_ids.fillna(method='bfill')
        cluster_ids = list(cluster_ids)

        # print cell_ids
        # print cluster_ids
        return cluster_ids

    @staticmethod
    def milen_cluster_larval_param_fit(cluster_id, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
        # Load Milen's best-fit larval params JSON file and get larval param fit for this cluster_id
        fn = base + "data/zambia/milen_clusters/milen_best_fits.json"
        f = open(fn, "r")
        larval_fits_dict = json.load(f)
        f.close()

        return_dict = larval_fits_dict[cluster_id]['fit']['params'].copy()
        del return_dict['drug_cov']
        del return_dict['itn_level']

        return return_dict

    @staticmethod
    def find_milen_larval_param_fit_for_grid_cells(cell_ids, fudge_milen_habitats=False):
        # Return the best-fit larval parameters for a set of grid cell IDs
        cluster_ids = ZambiaExperiment.find_milen_cluster_for_grid_cells(cell_ids)

        arab_params = np.zeros_like(cluster_ids)
        funest_params = np.zeros_like(cluster_ids)

        i = 0
        for cl_id in cluster_ids:
            # print "cl_id ",cl_id
            param_dict = ZambiaExperiment.milen_cluster_larval_param_fit(cl_id)
            arab_params[i] = param_dict['arabiensis_sc']
            funest_params[i] = param_dict['funestus_sc']

            if fudge_milen_habitats:
                [arab_params[i],funest_params[i]]= ZambiaExperiment.fudge_milen_habitats_by_hand(cl_id,
                                                                                 arab_params[i].astype(np.float),
                                                                                 funest_params[i].astype(np.float))

            i += 1

        arab_params = arab_params.astype(np.float)
        funest_params = funest_params.astype(np.float)
        return [arab_params, funest_params]

    @staticmethod
    def fudge_milen_habitats_by_hand(cl_id,arab_param, funest_param):
        # Modify larval habitat parameters by hand
        mod_factor_dict = {
            # Bbondo:
            "80201_1": 0.5,
            "80201_2": 0.5,
            "80201_5": 0.5,
            "80201_6": 0.5,
            "80201_8": 0.5,
            "80201_10": 0.5,
            "80201_12": 0.5,
            # Chabbobboma:
            "80202_2": 2.0,
            "80202_4": 2.0,
            "80202_7": 2.0,
            "80202_8": 2.0,
            "80202_9": 2.0,
            "80202_10": 2.0,
            "80203_6": 2.0,
            "80203_7": 2.0,
            # Chisanga:
            "80204_1": 1.5,
            "80204_4": 2.0,
            "80204_6": 1.5,
            "80204_7": 1.5,
            "80204_8": 2.0,
            "80204_9": 1.5,
            "80204_10": 1.5,
            # Chiyabi:
            "81102_1": 2.0,
            "81102_2": 2.0,
            "81102_3": 2.0,
            "81102_4": 2.0,
            "81102_5": 2.0,
            "81102_6": 2.0,
            "81102_9": 2.0,
            # Luumbo:
            "80208_1": 2.0,
            "80208_4": 2.0,
            "80208_5": 2.0,
            "80208_7": 2.0,
            "80208_8": 2.0,
            "80208_10": 2.0,
            # Munyumbwe:
            "80209_2": 1.5,
            "80209_5": 0.67,
            "80209_6": 2.0,
            "80209_7": 1.5,
            "80209_8": 0.67,
            # Nyanga Chamwe:
            "80210_3": 0.15,
            "80210_4": 0.15,
            "80210_5": 0.15,
            "80210_6": 0.15,
            "80210_7": 0.15,
            "80210_8": 0.15,
            "80210_9": 0.15,
            "80210_10": 0.15,
            "80210_11": 0.15,
            "80210_12": 0.15,
            # Sinafala:
            "80204_7": 2.0,
            "80211_1": 2.0,
            "80211_3": 2.0,
            "80211_4": 2.0,
            "80211_5": 1.5,
            "80211_6": 2.0,
            "80211_7": 2.0,
            "80211_8": 1.5,
            # Sinamalima:
            "81111_1": 2.0,
            "81111_2": 2.0,
            "81111_3": 2.0,
            "81111_4": 1.5,
            "81111_5": 2.0,
            "81111_6": 2.0,
            "81111_7": 1.5,
            "81111_8": 2.0,
        }
        try:
            return [arab_param * mod_factor_dict[cl_id],
                    funest_param * mod_factor_dict[cl_id]]
        except KeyError:
            return [arab_param, funest_param]

    @staticmethod
    def milen_cluster_climate_category(cluster_id, base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
        # Return the climate category for a given cluster ID
        lookup = pd.read_csv(base + "data/milen_clusters/cluster_climate_mapping.csv")
        return lookup[lookup['cluster'] == cluster_id]['climate.cat'].item()

    @staticmethod
    def closest_milen_immunity_overlay_filename(arab_param, funest_param, climate_category,
                                                base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
        # Find filename of immunity file which is closest match, within given climate category, to specified arabiensis/funestus larval habitat params
        def find_nearest_from_precomputed(x, precomputed_x):
            precomp_x_arr = np.array(precomputed_x)  # assume precomputed_x is a list
            closest = np.argmin(np.abs(x - precomp_x_arr))
            return precomputed_x[closest]

        # Get arab and funest precomputed values for specified climate category
        precomputed_arab_params, precomputed_funest_params = ZambiaExperiment.get_arab_funest_precomputed_values_for_climate_category(
            climate_category)

        nearest_arab_param = find_nearest_from_precomputed(arab_param, precomputed_arab_params)
        nearest_funest_param = find_nearest_from_precomputed(funest_param, precomputed_funest_params)

        # Return filename for these parameters:
        fn = base + "data/immunity/{}_1_node/{}_1_node_immune_init_p1_{}_p2_{}.json".format(climate_category,
                                                                                            climate_category,
                                                                                            nearest_funest_param,
                                                                                            nearest_arab_param)
        return fn

    @staticmethod
    def get_arab_funest_precomputed_values_for_climate_category(climate_category,
                                                                base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
        # Get all filenames in that climate folder, and extract the precomputed vector param values by parsing the filenames

        def represents_int(s):
            try:
                int(s)
                return True
            except ValueError:
                return False

        def get_arab_funest_from_filename(fn):
            # Extract the two parameters from the file name:
            break1 = fn.split('_p1_')
            break1 = break1[1]

            break2 = break1.split('_p2_')
            if represents_int(break2[0]):
                funest_param = int(break2[0])
            else:
                funest_param = np.float(break2[0])

            break2 = break2[1]
            break3 = break2.split('.json')
            if represents_int(break3[0]):
                arab_param = int(break3[0])
            else:
                arab_param = np.float(break3[0])

            return [arab_param, funest_param]

        import glob
        folder_name = base + "data/immunity/{}_1_node".format(climate_category)
        file_list = glob.glob(folder_name + "/*.json")

        arab_params = []  # np.zeros(len(file_list))
        funest_params = []  # np.zeros(len(file_list))

        # Loop over all filenames and extract vector parameters from the filenames
        # i = 0
        # print file_list
        for fn in file_list:
            a, f = get_arab_funest_from_filename(fn)
            # arab_params[i] = a
            # funest_params[i] = f
            arab_params.append(a)
            funest_params.append(f)

        # return [np.unique(arab_params),np.unique(funest_params)]
        unique_arab_params = list(set(arab_params))
        unique_funest_params = list(set(funest_params))

        unique_arab_params.sort()
        unique_funest_params.sort()
        return [unique_arab_params, unique_funest_params]

    @staticmethod
    def closest_milen_immunity_overlay_filenames_for_grid_cells(cell_ids):
        cluster_ids = ZambiaExperiment.find_milen_cluster_for_grid_cells(cell_ids)

        fn_list = []
        for cluster_id in cluster_ids:
            larval_param_dict = ZambiaExperiment.milen_cluster_larval_param_fit(cluster_id)
            climate_category = ZambiaExperiment.milen_cluster_climate_category(cluster_id)

            fn = ZambiaExperiment.closest_milen_immunity_overlay_filename(larval_param_dict['arabiensis_sc'],
                                                         larval_param_dict['funestus_sc'],
                                                         climate_category)

            fn_list.append(fn)

        return fn_list