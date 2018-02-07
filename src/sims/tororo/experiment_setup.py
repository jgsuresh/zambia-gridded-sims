"""
Run multi-node simulation in Lake Kariba region.
"""

import os
import json
import pandas as pd
import numpy as np

# import matplotlib
# matplotlib.use('')

from dtk.tools.climate.ClimateGenerator import ClimateGenerator
from dtk.tools.migration.MigrationGenerator import MigrationGenerator
from dtk.tools.spatialworkflow.DemographicsGenerator import DemographicsGenerator
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.vector.species import set_species_param
from dtk.interventions.itn_age_season import add_ITN_age_season
from dtk.interventions.health_seeking import add_health_seeking
from dtk.interventions.irs import add_IRS
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import ModBuilder, ModFn
from simtools.SetupParser import SetupParser
from malaria.interventions.malaria_drug_campaigns import add_drug_campaign
# from dtk.utils.reports.MalariaReport import add_summary_report
from malaria.reports.MalariaReport import add_summary_report

from relative_time import *
from gridded_sim_general import *
from grid_ids_to_nodes import generate_lookup
from gen_migr_json import gen_gravity_links_json, save_link_rates_to_txt
from gridded_sim_general import *



class COMPS_Experiment:

    def __init__(self,
                 base,
                 exp_name,
                 catch='all',
                 grid_pop_csv_file=None,
                 imm_1node_fp=None,
                 larval_params_mode='uniform',
                 immunity_mode='uniform',
                 migration_on=True,
                 gravity_migr_params=None,
                 rcd_people_num=5,
                 start_year=2001,
                 sim_length_years=19,
                 num_cores=12,
                 healthseek_fn=None,
                 itn_fn=None,
                 irs_fn=None,
                 msat_fn=None,
                 mda_fn=None,
                 stepd_fn=None,
                 fudge_milen_habitats=False):

        """
        :param base:
        :param exp_name:
        :param catch:
        :param grid_pop_csv_file:
        :param imm_1node_fp: only necessary if immunity_mode is set to "uniform"
        :param larval_params_mode:
        :param immunity_mode: "naive" if no immunity initialization.  "uniform" to use immunity file from single node for all nodes.  "milen" to draw immunity initializations from Milen's immunity overlay grids
        :param migration_on:
        :param gravity_migr_params:
        :param rcd_people_num:
        :param start_year:
        :param sim_length_years:
        :param num_cores:
        :param healthseek_fn:
        :param itn_fn:
        :param irs_fn:
        :param msat_fn:
        :param mda_fn:
        :param stepd_fn:
        :param fudge_milen_habitats:
        """

        self.base = base
        self.exp_name = exp_name
        self.exp_base = base + 'data/COMPS_experiments/{}/'.format(exp_name)
        self.catch = catch
        self.grid_pop_csv_file = grid_pop_csv_file
        self.imm_1node_fp = imm_1node_fp
        self.larval_params_mode = larval_params_mode
        self.immunity_mode = immunity_mode

        self.migration_on = migration_on
        self.gravity_migr_params = gravity_migr_params
        self.rcd_people_num = rcd_people_num

        self.start_year = start_year
        self.sim_length_years = sim_length_years

        self.num_cores = num_cores

        self.healthseek_fn = healthseek_fn
        self.itn_fn = itn_fn
        self.irs_fn = irs_fn
        self.msat_fn = msat_fn
        self.mda_fn = mda_fn
        self.stepd_fn = stepd_fn

        self.fudge_milen_habitats = fudge_milen_habitats

        # Ensure directories exist:
        self.ensure_filesystem()

        self.cb = self.build_cb()
        self.multinode_setup() # Assume no single-node for now
        self.basic_sim_setup()


    #################################################################################################
    # BASIC SETUP (CONFIG BUILDER + DEMOGRAPHICS FILE)

    def ensure_filesystem(self):
        def ensure_dir(file_path):
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

        ensure_dir(self.exp_base)
        ensure_dir(self.exp_base + 'Demographics/')
        ensure_dir(self.exp_base + 'Immunity/')
        ensure_dir(self.exp_base + 'Climate/')
        ensure_dir(self.exp_base + 'Migration/')
        ensure_dir(self.exp_base + 'Logs/')

    def build_cb(self):
        location = 'HPC'
        SetupParser.default_block = location

        cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')
        cb.set_experiment_executable(self.base + 'bin/Eradication.exe')
        cb.set_input_files_root(self.exp_base)

        # Tell config builder where to find dlls for specified Bamboo build of executable
        cb.set_dll_root(self.base + 'bin/')

        return cb

    def basic_sim_setup(self, record_events=True):

        # Miscellaneous:
        self.cb.set_param("Enable_Demographics_Other", 1)
        self.cb.set_param("Enable_Demographics_Builtin", 0)
        self.cb.set_param("Valid_Intervention_States", [])
        self.cb.set_param("New_Diagnostic_Sensitivity", 0.025) # 40/uL

        # Human population properties:
        self.cb.update_params({
            'Birth_Rate_Dependence': 'FIXED_BIRTH_RATE',  # Match demographics file for constant population size (with exponential age distribution)
            'Enable_Nondisease_Mortality': 1,
        })


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


        # Immunity:
        if self.immunity_mode == "naive":
            self.cb.update_params({
                "Enable_Immunity_Initialization_Distribution": 0
            })
        else:
            self.cb.update_params({
                    "Enable_Immunity_Initialization_Distribution": 1,
                    "Immunity_Initialization_Distribution_Type": "DISTRIBUTION_COMPLEX",
            })


        # Event recording:
        intervene_events_list = ["Bednet_Got_New_One","Bednet_Using","Bednet_Discarded"]
        migration_events_list = ["Immigrating", "Emigrating"]
        # case_events_list = ["NewClinicalCase",

        full_events_list = intervene_events_list # migration events too verbose

        self.cb.update_params({
            "Report_Event_Recorder": 1,
            "Report_Event_Recorder_Ignore_Events_In_List": 0,
            "Listed_Events": full_events_list,
            "Report_Event_Recorder_Events": full_events_list
        })


        # Spatial reporting:
        self.cb.update_params({
            'Enable_Spatial_Output': 1,  # turn on spatial reporting
            'Spatial_Output_Channels': ['Infectious_Vectors', 'Adult_Vectors', 'New_Infections', 'Population',
                                        'Prevalence',
                                        'New_Diagnostic_Prevalence', 'Daily_EIR', 'New_Clinical_Cases',
                                        'Human_Infectious_Reservoir', 'Daily_Bites_Per_Human',
                                        'Land_Temperature',
                                        'Relative_Humidity', 'Rainfall', 'Air_Temperature']
        })

        add_summary_report(self.cb)

        # Basic health-seeking
        self.implement_dummy_healthseeking()

    def implement_dummy_healthseeking(self):
        # Implement dummy health-seeking behavior for one day at beginning of simulation, to avoid DrugStatus error.
        start_date = "{}-01-01".format(self.start_year)  # Day 1 of simulation
        date_format = "%Y-%m-%d"
        sim_duration = 365 * self.sim_length_years  # length in days
        self.cb.params['Simulation_Duration'] = sim_duration

        start_date_refmt = convert_to_day(start_date, start_date, date_format)

        add_health_seeking(self.cb,
                           start_day=float(start_date_refmt),
                           targets=[{'trigger': 'NewClinicalCase',
                                     'coverage': 0.1,
                                     'agemin': 0,
                                     'agemax': 5,
                                     'seek': 1, 'rate': 0.3},
                                    {'trigger': 'NewClinicalCase',
                                     'coverage': 0.1,
                                     'agemin': 5,
                                     'agemax': 100,
                                     'seek': 1,
                                     'rate': 0.3},
                                    {'trigger': 'NewSevereCase',
                                     'coverage': 0.1,
                                     'agemin': 0,
                                     'agemax': 5,
                                     'seek': 1, 'rate': 0.3},
                                    {'trigger': 'NewSevereCase',
                                     'coverage': 0.1,
                                     'agemin': 5,
                                     'agemax': 100,
                                     'seek': 1,
                                     'rate': 0.3}],
                           drug=['Artemether', 'Lumefantrine'],
                           dosing='FullTreatmentNewDetectionTech',
                           nodes={"class": "NodeSetAll"},
                           duration=float(1))

    def implement_baseline_healthseeking(self):
        # Implement basic health-seeking behavior for all individuals in simulation
        start_date = "{}-01-01".format(self.start_year)  # Day 1 of simulation
        date_format = "%Y-%m-%d"
        sim_duration = 365 * self.sim_length_years  # length in days
        self.cb.params['Simulation_Duration'] = sim_duration

        # Prevent DTK from spitting out too many messages
        self.cb.params['logLevel_JsonConfigurable'] = "WARNING"
        self.cb.params['Disable_IP_Whitelist'] = 1

        # Event information files

        if self.healthseek_fn==None:
            healthseek_event_file = self.base + 'data/interventions/chiyabi/gridded-uniform/grid_chiyabi_hfca_healthseek_events.csv'
        else:
            healthseek_event_file = self.healthseek_fn
        healthseek_events = pd.read_csv(healthseek_event_file)

        # Compute simulation days relative to start date or use default in file
        healthseek_events['simday'] = [convert_to_day(x, start_date, date_format) for x in
                                       healthseek_events.fulldate]

        # Get grid cell to node ID lookup table:
        nodeid_lookup,pop_lookup = generate_lookup(self.demo_fp_full)

        # Restrict to catchment of interest
        catch_cells = find_cells_for_this_catchment(self.catch)
        healthseek_events = healthseek_events[np.in1d(healthseek_events['grid_cell'],catch_cells)]
        healthseek_events.reset_index(inplace=True)

        for hs in range(len(healthseek_events)):
            node_list = [nodeid_lookup[healthseek_events['grid_cell'][hs]]]

            add_health_seeking(self.cb,
                               start_day=float(healthseek_events['simday'][hs]),
                               targets=[{'trigger': 'NewClinicalCase',
                                         'coverage': float(healthseek_events['cov_newclin_youth'][hs]), 'agemin': 0,
                                         'agemax': 5,
                                         'seek': 1, 'rate': 0.3},
                                        {'trigger': 'NewClinicalCase',
                                         'coverage': float(healthseek_events['cov_newclin_adult'][hs]), 'agemin': 5,
                                         'agemax': 100,
                                         'seek': 1,
                                         'rate': 0.3},
                                        {'trigger': 'NewSevereCase',
                                         'coverage': float(healthseek_events['cov_severe_youth'][hs]), 'agemin': 0,
                                         'agemax': 5,
                                         'seek': 1, 'rate': 0.3},
                                        {'trigger': 'NewSevereCase',
                                         'coverage': float(healthseek_events['cov_severe_adult'][hs]), 'agemin': 5,
                                         'agemax': 100,
                                         'seek': 1,
                                         'rate': 0.3}],
                               drug=['Artemether', 'Lumefantrine'],
                               dosing='FullTreatmentNewDetectionTech',
                               nodes={"class": "NodeSetNodeList", "Node_List": node_list},
                               duration=float(healthseek_events['duration'][hs]))

    def multinode_setup(self):
        self.demo_fp_from_input = "Demographics/demo.json"
        self.demo_fp_full = os.path.join(self.exp_base,self.demo_fp_from_input)
        self.cb.update_params({'Demographics_Filenames': [self.demo_fp_from_input]})

        if self.immunity_mode != "naive":
            self.immun_fp_from_input = "Immunity/immun.json"
            self.immun_fp_full = os.path.join(self.exp_base, self.immun_fp_from_input)
            self.cb.update_params({'Demographics_Filenames': [self.demo_fp_from_input, self.immun_fp_from_input]})

        self.cb.set_param("Num_Cores", self.num_cores)

        #######################################################################################################
        # CLIMATE-RELATED PARAMETERS:
        #######################################################################################################
        self.cb.update_params({
            'Air_Temperature_Filename': "Climate/Zambia_30arcsec_air_temperature_daily.bin",
            'Land_Temperature_Filename': "Climate/Zambia_30arcsec_air_temperature_daily.bin",
            'Rainfall_Filename': "Climate/Zambia_30arcsec_rainfall_daily.bin",
            'Relative_Humidity_Filename': "Climate/Zambia_30arcsec_relative_humidity_daily.bin"
        })

        #######################################################################################################
        # MIGRATION-RELATED PARAMETERS:
        #######################################################################################################
        if self.migration_on:
            self.cb.update_params({
                                    'Migration_Model': 'FIXED_RATE_MIGRATION',
                                    'Local_Migration_Filename': 'Migration/local_migration.bin', # note that underscore prior 'migration.bin' is required for legacy reasons that need to be refactored...
                                    'Enable_Local_Migration':1,
                                    'Migration_Pattern': 'SINGLE_ROUND_TRIPS', # human migration
                                    'Local_Migration_Roundtrip_Duration': 2, # mean of exponential days-at-destination distribution
                                    'Local_Migration_Roundtrip_Probability': 0.95, # fraction that return
                                    'x_Local_Migration': 4 # amplitude of migration
            })
        elif not self.migration_on:
            self.cb.update_params({'Migration_Model': 'NO_MIGRATION'})  #'NO_MIGRATION' is actually default for MALARIA_SIM, but might as well make sure it's off

    def gen_demo_file(self,input_csv,larval_params=None):
        dg = DemographicsGenerator.from_file(self.cb, input_csv, catch=self.catch)
        demo_dict = dg.generate_demographics()

        # Add catchment name to demographics file metadata:
        demo_dict["Metadata"]["Catchment"] = self.catch

        # Add larval habitat parameters to demographics file:
        if self.larval_params_mode != "calibrate":
            demo_dict = self.add_larval_habitats_to_demo(demo_dict, larval_params=larval_params)

        if larval_params:
            temp_h = larval_params['temp_h']
            linear_h = larval_params['linear_h']
            demo_fp = self.exp_base + "Demographics/demo_temp{}_linear{}.json".format(int(temp_h),int(linear_h))
        else:
            demo_fp = self.exp_base + "Demographics/demo.json"

        demo_f = open(demo_fp, 'w+')
        json.dump(demo_dict, demo_f, indent=4)
        demo_f.close()


    #################################################################################################
    # LARVAL PARAMETER-RELATED FUNCTIONS:
    def add_larval_habitats_to_demo(self, demo_dict, nodeset='all', larval_params=None):
        # Add larval habitat parameters (some of which scale with population) to demographics file

        def add_larval_habitat_to_node(node_item, const_h,temp_h,water_h,linear_h):
            calib_single_node_pop = 1000  # for Zambia

            # This is now done in the demographics generator itself:
            # birth_rate = (float(node_item['NodeAttributes']['InitialPopulation']) / (1000 + 0.0)) * 0.12329
            # node_item['NodeAttributes']['BirthRate'] = birth_rate

            pop_multiplier = float(node_item['NodeAttributes']['InitialPopulation']) / (calib_single_node_pop + 0.0)

            temp_multiplier = temp_h * pop_multiplier
            linear_multiplier = linear_h * pop_multiplier
            const_multiplier = const_h # NOTE: No pop multiplier
            water_multiplier = water_h * pop_multiplier

            node_item['NodeAttributes']['LarvalHabitatMultiplier'] = {
                "CONSTANT": const_multiplier,
                "TEMPORARY_RAINFALL": temp_multiplier,
                "WATER_VEGETATION": water_multiplier,
                "LINEAR_SPLINE": linear_multiplier
            }


        if self.larval_params_mode=='milen':
            # get grid cells from pop csv file:
            # for those grid cells, get corresponding arab/funest params
            # loop over nodes [order will correspond, by construction, to pop csv ordering]
            # give each node the corresponding larval params

            # Load pop csv file to get grid cell numbers:
            # pop_df = pd.read_csv(self.grid_pop_csv_file)
            # grid_cells = np.array(pop_df['node_label'])
            grid_cells = find_cells_for_this_catchment(self.catch)

            # From those grid cells, and the Milen-clusters they correspond to, get best-fit larval habitat parameters
            arab_params,funest_params = find_milen_larval_param_fit_for_grid_cells(grid_cells,fudge_milen_habitats=self.fudge_milen_habitats)

            # Loop over nodes in demographics file (which will, by construction, correspond to the grid pop csv ordering)
            i = 0
            for node_item in demo_dict['Nodes']:
                # if larval_params:
                #     const_h = larval_params['const_h']
                #     temp_h = larval_params['temp_h']
                #     water_h = larval_params['water_h']
                #     linear_h = larval_params['linear_h']
                # else:
                const_h = 1.
                temp_h = arab_params[i]
                water_h = 1.
                linear_h = funest_params[i]

                add_larval_habitat_to_node(node_item,const_h,temp_h,water_h,linear_h)

                i += 1

        elif self.larval_params_mode=='uniform':
            for node_item in demo_dict['Nodes']:

                if larval_params:
                    const_h = larval_params['const_h']
                    temp_h = larval_params['temp_h']
                    water_h = larval_params['water_h']
                    linear_h = larval_params['linear_h']
                else:
                    const_h = 1.
                    temp_h = 122.
                    water_h = 1.
                    linear_h = 97.

                add_larval_habitat_to_node(node_item,const_h,temp_h,water_h,linear_h)

        return demo_dict

    def larval_param_sweeper(self,cb,temp_h,linear_h):
        larval_params = {
            "const_h": 1.,
            "water_h": 1.,
            "temp_h": np.float(temp_h),
            "linear_h": np.float(linear_h)
        }

        # Will need new demographics file that incorporates these larval parameters.
        # demographics filename example: Demographics/MultiNode/demo_temp50_linear120.json
        new_demo_fp_from_input = self.demo_fp_from_input[:-5] + "_temp{}_linear{}.json".format(temp_h, linear_h)

        # Check if demographics file for these parameters already exists.  If not, create it
        if not os.path.isfile(new_demo_fp_from_input):
            self.gen_demo_file(self.grid_pop_csv_file,larval_params=larval_params)

        # Then pass this demographics file to the config_builder.
        if self.immunity_mode == "naive":
            self.cb.update_params({'Demographics_Filenames': [new_demo_fp_from_input]})
        else:
            self.cb.update_params({'Demographics_Filenames': [new_demo_fp_from_input, self.immun_fp_from_input]})


        return larval_params

    #################################################################################################
    # MIGRATION-RELATED FUNCTIONS:
    def gen_migration_files(self):
        migr_json_fp = self.exp_base + "Migration/grav_migr_rates.json"

        migr_dict = gen_gravity_links_json(self.demo_fp_full, self.gravity_migr_params, outf=migr_json_fp)
        rates_txt_fp = self.exp_base + "Migration/grav_migr_rates.txt"

        save_link_rates_to_txt(rates_txt_fp, migr_dict)

        # Generate migration binary:
        migration_filename = self.cb.get_param('Local_Migration_Filename')
        print("migration_filename: ",migration_filename)
        MigrationGenerator.link_rates_txt_2_bin(rates_txt_fp,
                                                self.exp_base+migration_filename)

        # Generate migration header:
        MigrationGenerator.save_migration_header(self.demo_fp_full,
                                                 outfilename=self.exp_base +'Migration/local_migration.bin.json'
                                                 )

    def vector_migration_sweeper(self, vector_migration_on):
        if vector_migration_on:
            self.cb.update_params({
                'Vector_Migration_Modifier_Equation': 'LINEAR',
                'Vector_Sampling_Type': 'SAMPLE_IND_VECTORS', # individual vector model (required for vector migration)
                'Mosquito_Weight': 10,
                'Enable_Vector_Migration': 1, # mosquito migration
                'Enable_Vector_Migration_Local': 1, # migration rate hard-coded in NodeVector::processEmigratingVectors() such that 50% total leave a 1km x 1km square per day (evenly distributed among the eight adjacent grid cells).
                'Vector_Migration_Base_Rate': 0.15, # default is 0.5
                'x_Vector_Migration_Local': 1
            })
        else:
            self.cb.update_params({
                'Enable_Vector_Migration': 0,  # mosquito migration
                'Enable_Vector_Migration_Local': 0
            # migration rate hard-coded in NodeVector::processEmigratingVectors() such that 50% total leave a 1km x 1km square per day (evenly distributed among the eight adjacent grid cells).
            })
        return {"vec_migr": vector_migration_on}

    #################################################################################################
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

    def implement_interventions(self, cb, include_itn, include_irs, include_msat, include_mda, include_stepd):
        start_date = "{}-01-01".format(self.start_year)  # Day 1 of simulation
        date_format = "%Y-%m-%d"
        sim_duration = 365 * self.sim_length_years  # length in days
        self.cb.params['Simulation_Duration'] = sim_duration

        # Prevent DTK from spitting out too many messages
        self.cb.params['logLevel_JsonConfigurable'] = "WARNING"
        self.cb.params['Disable_IP_Whitelist'] = 1

        # Event information files
        if self.itn_fn == None:
            itn_event_file = self.base + 'data/interventions/chiyabi/gridded-uniform/grid_chiyabi_hfca_itn_events.csv'
        else:
            itn_event_file = self.itn_fn

        if self.irs_fn == None:
            irs_event_file = self.base + 'data/interventions/chiyabi/gridded-uniform/grid_chiyabi_hfca_irs_events.csv'
        else:
            irs_event_file = self.irs_fn

        if self.msat_fn == None:
            msat_event_file = self.base + 'data/interventions/chiyabi/gridded-uniform/grid_chiyabi_hfca_msat_events.csv'
        else:
            msat_event_file = self.msat_fn

        if self.mda_fn == None:
            mda_event_file = self.base + 'data/interventions/chiyabi/gridded-uniform/grid_chiyabi_hfca_mda_events.csv'
        else:
            mda_event_file = self.mda_fn

        if self.stepd_fn == None:
            stepd_event_file = self.base + 'data/interventions/chiyabi/gridded-uniform/grid_chiyabi_hfca_stepd_events.csv'
        else:
            stepd_event_file = self.stepd_fn

        # Import event info
        itn_events = pd.read_csv(itn_event_file)
        irs_events = pd.read_csv(irs_event_file)
        msat_events = pd.read_csv(msat_event_file)
        mda_events = pd.read_csv(mda_event_file)
        stepd_events = pd.read_csv(stepd_event_file)


        # Compute simulation days relative to start date or use default in file
        itn_events['simday'] = [convert_to_day(x, start_date, date_format) for x in itn_events.fulldate]
        irs_events['simday'] = [convert_to_day(x, start_date, date_format) for x in irs_events.fulldate]
        msat_events['simday'] = [convert_to_day(x, start_date, date_format) for x in msat_events.fulldate]
        mda_events['simday'] = [convert_to_day(x, start_date, date_format) for x in mda_events.fulldate]
        stepd_events['simday'] = [convert_to_day(x, start_date, date_format) for x in stepd_events.fulldate]

        # Restrict to catchment of interest
        catch_cells = find_cells_for_this_catchment(self.catch)
        itn_events = itn_events[np.in1d(itn_events['grid_cell'], catch_cells)]
        irs_events = irs_events[np.in1d(irs_events['grid_cell'], catch_cells)]
        msat_events = msat_events[np.in1d(msat_events['grid_cell'], catch_cells)]
        mda_events = mda_events[np.in1d(mda_events['grid_cell'], catch_cells)]
        stepd_events = stepd_events[np.in1d(stepd_events['grid_cell'], catch_cells)]

        itn_events.reset_index(inplace=True)
        irs_events.reset_index(inplace=True)
        msat_events.reset_index(inplace=True)
        mda_events.reset_index(inplace=True)
        stepd_events.reset_index(inplace=True)

        # Get grid cell to node ID lookup table:
        from grid_ids_to_nodes import generate_lookup
        nodeid_lookup,pop_lookup = generate_lookup(self.demo_fp_full)


        if include_itn:
            for itn in range(len(itn_events)):

                # Add non-birth nets
                add_ITN_age_season(cb, start=float(itn_events['simday'][itn]),
                                   age_dep={'youth_cov': float(itn_events['age_cov'][itn]), 'youth_min_age': 5,
                                            'youth_max_age': 20},
                                   coverage_all=float(itn_events['cov_all'][itn]),
                                   as_birth=False,
                                   seasonal_dep={'min_cov': float(itn_events['min_season_cov'][itn]), 'max_day': 60},
                                   discard={'halflife1': 260, 'halflife2': 2106,
                                            'fraction1': float(itn_events['fast_fraction'][itn])},
                                   nodeIDs=[nodeid_lookup[itn_events['grid_cell'][itn]]])
                # Add birth nets
                if itn < (len(itn_events)-1) and (itn_events['grid_cell'][itn+1] == itn_events['grid_cell'][itn]):
                    birth_duration = float(itn_events['simday'][itn + 1] - itn_events['simday'][itn] - 1)
                else:
                    birth_duration = -1

                # This code only works for non-pixel-specific file format:
                # if itn < (len(itn_events) - 1):
                #     birth_duration = float(itn_events['simday'][itn + 1] - itn_events['simday'][itn] - 1)
                # else:
                #     birth_duration = -1

                add_ITN_age_season(cb, start=float(itn_events['simday'][itn]),
                                   age_dep={'youth_cov': float(itn_events['age_cov'][itn]), 'youth_min_age': 5,
                                            'youth_max_age': 20},
                                   coverage_all=float(itn_events['cov_all'][itn]),
                                   as_birth=True,
                                   seasonal_dep={'min_cov': float(itn_events['min_season_cov'][itn]), 'max_day': 60},
                                   discard={'halflife1': 260, 'halflife2': 2106,
                                            'fraction1': float(itn_events['fast_fraction'][itn])},
                                   duration=birth_duration,
                                   nodeIDs=[nodeid_lookup[itn_events['grid_cell'][itn]]])

        if include_irs:
            for irs in range(len(irs_events)):
                add_IRS(cb, start=float(irs_events['simday'][irs]),
                        coverage_by_ages=[{'coverage': float(irs_events['cov_all'][irs])}],
                        initial_killing=float(irs_events['killing'][irs]), duration=float(irs_events['duration'][irs]),
                        nodeIDs=[nodeid_lookup[irs_events['grid_cell'][irs]]])


        if include_msat:
            for msat in range(len(msat_events)):
                    add_drug_campaign(cb, campaign_type='MSAT', drug_code='AL',
                                      start_days=[float(msat_events['simday'][msat])],
                                      coverage=msat_events['cov_all'][msat], repetitions=1, interval=60,
                                      dosing='SingleDose',
                                      nodes=[nodeid_lookup[msat_events['grid_cell'][msat]]])

        if include_mda:
            for mda in range(len(mda_events)):
                    add_drug_campaign(cb, campaign_type='MDA', drug_code='DP',
                                      start_days=[float(mda_events['simday'][mda])],
                                      coverage=float(mda_events['cov_all'][mda]), repetitions=1, interval=60,
                                      dosing='SingleDose',
                                      nodes=[nodeid_lookup[mda_events['grid_cell'][mda]]])

        if include_stepd:
            for sd in range(len(stepd_events)):
                cov = np.min([1.,float(self.rcd_people_num) / float(pop_lookup[stepd_events['grid_cell'][sd]])])
                add_drug_campaign(cb, campaign_type='rfMDA', drug_code='AL',
                                  start_days=[float(stepd_events['simday'][sd])],
                                  # coverage=float(stepd_events['coverage'][sd]),
                                  # coverage=float(self.rcd_people_num)/float(self.pop_start),
                                  coverage=cov,
                                  interval=float(stepd_events['interval'][sd]),
                                  nodes=[nodeid_lookup[stepd_events['grid_cell'][sd]]])

        return {"ITNs": include_itn,
                "IRS": include_irs,
                "MSAT": include_msat,
                "MDA": include_mda,
                "StepD": include_stepd}

    def file_setup(self,generate_immunity_file=True,generate_demographics_file=True,generate_climate_files=True,generate_migration_files=True):
        # self.multinode_setup()

        if generate_demographics_file:
            print("Generating demographics file...")
            self.gen_demo_file(self.grid_pop_csv_file)

        if generate_immunity_file:
            print("Generating immunity files...")
            if self.immunity_mode == "uniform":
                self.get_immunity_file_from_single_node()
            elif self.immunity_mode == "milen":
                self.gen_immunity_file_from_milen_clusters()

        if generate_migration_files:
            print("Generating migration files...")
            self.gen_migration_files()

        if generate_climate_files:
            print("Generating climate files...")
            SetupParser.init()
            cg = ClimateGenerator(self.demo_fp_full,
                                  self.exp_base + 'Logs/climate_wo.json',
                                  self.exp_base + 'Climate/',
                                  start_year = str(self.start_year),
                                  num_years = str(np.min([2016 - self.start_year, self.sim_length_years]))
                                  )

            cg.generate_climate_files()



    def submit_experiment(self,num_seeds=1,intervention_sweep=False,larval_sweep=False,migration_sweep=False,vector_migration_sweep=False,
                          simple_intervention_sweep=True,custom_name=None):

        # Implement the actual (not dummy) baseline healthseeking
        self.implement_baseline_healthseeking()


        modlists = []

        if num_seeds > 1:
            new_modlist = [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed) for seed in range(num_seeds)]
            modlists.append(new_modlist)

        if larval_sweep:
            new_modlist = [ModFn(self.larval_param_sweeper, temp_h, linear_h)
                           for temp_h in [61,122,244]
                           for linear_h in [48, 97, 194]]
            modlists.append(new_modlist)

        if migration_sweep:
            new_modlist = [ModFn(DTKConfigBuilder.set_param, 'x_Local_Migration', x) for x in [0.5,1,2,5]]
            modlists.append(new_modlist)

        if vector_migration_sweep:
            new_modlist = [ModFn(self.vector_migration_sweeper, vector_migration_on) for vector_migration_on in [True, False]]
            modlists.append(new_modlist)

        if simple_intervention_sweep:
            new_modlist = [
                ModFn(self.implement_interventions, True, False, False, False, False),
                ModFn(self.implement_interventions, False, True, False, False, False),
                ModFn(self.implement_interventions, False, False, True, False, False),
                ModFn(self.implement_interventions, False, False, False, True, False),
                ModFn(self.implement_interventions, False, False, False, False, True),
                ModFn(self.implement_interventions, True, True, True, True, True)
            ]
            modlists.append(new_modlist)
        else:
            new_modlist = [ModFn(self.implement_interventions, True, True, True, True, True)]
            modlists.append(new_modlist)



        # if intervention_sweep:
        #     # Interventions to turn on or off
        #     include_itn_list = [True, False]
        #     include_irs_list = [True, False]
        #     include_mda_list = [True, False]
        #     include_msat_list = [True, False]
        #     include_stepd_list = [True, False]
        #
        #     new_modlist = [
        #         ModFn(self.implement_interventions, use_itn, use_irs, use_msat, use_mda, use_stepd)
        #         for use_itn in include_itn_list
        #         for use_irs in include_irs_list
        #         for use_mda in include_mda_list
        #         for use_msat in include_msat_list
        #         for use_stepd in include_stepd_list
        #     ]
        #
        # else:
        #     new_modlist = [ModFn(self.implement_interventions, True, True, True, True, True)]
        # modlists.append(new_modlist)


        builder = ModBuilder.from_combos(*modlists)

        run_name = self.exp_name
        if custom_name:
            run_name = custom_name


        # SetupParser.init()
        # SetupParser.set("HPC","priority","Normal")
        exp_manager = ExperimentManagerFactory.init()
        exp_manager.run_simulations(config_builder=self.cb, exp_name=run_name, exp_builder=builder)


    ################################################################################################################
    def return_cb_for_calibration(self):
        self.implement_baseline_healthseeking()
        self.implement_interventions(self.cb,True, True, True, True, True)
        return self.cb