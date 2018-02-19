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

    def gen_immunity_file_from_milen_clusters(self):
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
            immun_fn_list = closest_milen_immunity_overlay_filenames_for_grid_cells(cell_ids)

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

    def return_cb_for_calibration(self):
        self.implement_baseline_healthseeking()
        self.implement_interventions(self.cb,True, True, True, True, True)
        return self.cb
