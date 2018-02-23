"""
Run multi-node simulation in Lake Kariba region.
"""
import csv
import os
import json
import pandas as pd
import numpy as np

import matplotlib

from dtk.tools.demographics.Node import Node, nodeid_from_lat_lon

matplotlib.use('Agg')

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
from malaria.reports.MalariaReport import add_summary_report

from relative_time import *
from gridded_sim_general import *
from grid_ids_to_nodes import generate_lookup
from gen_migr_json import gen_gravity_links_json, save_link_rates_to_txt

# from zambia_experiments import ZambiaExperiment


class GriddedConfigBuilder:

    def __init__(self,
                 base,
                 exp_name,
                 desired_cells,
                 region="Zambia",
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

        self.base = base
        self.exp_name = exp_name
        self.exp_base = base + 'data/COMPS_experiments/{}/'.format(exp_name)
        self.desired_cells = desired_cells
        self.region = region
        self.start_year = start_year
        self.sim_length_years = sim_length_years
        self.immunity_mode = immunity_mode
        self.num_cores = num_cores
        self.parser_location = parser_location

        self.healthseek_fn = healthseek_fn
        self.itn_fn = itn_fn
        self.irs_fn = irs_fn
        self.msat_fn = msat_fn
        self.mda_fn = mda_fn
        self.stepd_fn = stepd_fn


        self.demo_fp_from_input = "Demographics/demo.json"
        self.demo_fp_full = os.path.join(self.exp_base,self.demo_fp_from_input)

        self.ensure_filesystem() # Need this to build config-builder
        # Build config-builder
        self.build_cb()


    def ensure_filesystem(self):
        def ensure_dir(file_path):
            directory = os.path.dirname(file_path)
            if not os.path.exists(directory):
                os.makedirs(directory)

        ensure_dir(self.exp_base)
        ensure_dir(self.exp_base + 'Demographics/')
        # if self.immunity_mode != "naive":
        ensure_dir(self.exp_base + 'Immunity/')
        ensure_dir(self.exp_base + 'Climate/')
        ensure_dir(self.exp_base + 'Migration/')
        ensure_dir(self.exp_base + 'Logs/')

    #################################################################################################
    # CONFIG-BUILDER SETUP

    def build_cb(self):
        self.cb = self.basic_cb()
        self.spatial_cb_setup()
        self.add_reporting_to_cb()
        self.implement_dummy_healthseeking()

    def basic_cb(self):
        SetupParser.default_block = self.parser_location

        cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')
        cb.set_experiment_executable(self.base + 'bin/Eradication.exe')
        cb.set_input_files_root(self.exp_base)

        # Tell config builder where to find dlls for specified Bamboo build of executable
        cb.set_dll_root(self.base + 'bin/')

        cb.set_param("Num_Cores", self.num_cores)
        return cb

    def spatial_cb_setup(self, migration_on=True):
        self.cb.update_params({'Demographics_Filenames': [self.demo_fp_from_input]})

        if self.immunity_mode != "naive":
            self.immun_fp_from_input = "Immunity/immun.json"
            self.immun_fp_full = os.path.join(self.exp_base, self.immun_fp_from_input)
            self.cb.update_params({'Demographics_Filenames': [self.demo_fp_from_input, self.immun_fp_from_input]})


        # CLIMATE
        self.cb.update_params({
            'Air_Temperature_Filename': "Climate/{}_30arcsec_air_temperature_daily.bin".format(self.region),
            'Land_Temperature_Filename': "Climate/{}_30arcsec_air_temperature_daily.bin".format(self.region),
            'Rainfall_Filename': "Climate/{}_30arcsec_rainfall_daily.bin".format(self.region),
            'Relative_Humidity_Filename': "Climate/{}_30arcsec_relative_humidity_daily.bin".format(self.region)
        })
        #######################################################################################################
        # MIGRATION-RELATED PARAMETERS:
        #######################################################################################################
        if migration_on:
            self.cb.update_params({
                                'Migration_Model': 'FIXED_RATE_MIGRATION',
                                'Local_Migration_Filename': 'Migration/local_migration.bin', # note that underscore prior 'migration.bin' is required for legacy reasons that need to be refactored...
                                'Enable_Local_Migration':1,
                                'Migration_Pattern': 'SINGLE_ROUND_TRIPS', # human migration
                                'Local_Migration_Roundtrip_Duration': 2, # mean of exponential days-at-destination distribution
                                'Local_Migration_Roundtrip_Probability': 0.95, # fraction that return
            })
        else:
            self.cb.update_params({'Migration_Model': 'NO_MIGRATION'})  #'NO_MIGRATION' is actually default for MALARIA_SIM, but might as well make sure it's off

    def add_reporting_to_cb(self, record_events=False):

        # Miscellaneous:
        self.cb.set_param("Enable_Demographics_Other", 1)
        self.cb.set_param("Enable_Demographics_Builtin", 0)
        self.cb.set_param("Valid_Intervention_States", [])
        # self.cb.set_param("New_Diagnostic_Sensitivity", 0.025) # 40/uL
        self.cb.set_param("Report_Detection_Threshold_True_Parasite_Density", 40.0)

        # Human population properties:
        self.cb.update_params({
            'Birth_Rate_Dependence': 'FIXED_BIRTH_RATE',  # Match demographics file for constant population size (with exponential age distribution)
            'Enable_Nondisease_Mortality': 1,
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

        if record_events:
            self.cb.update_params({
                "Report_Event_Recorder": 1,
                "Report_Event_Recorder_Ignore_Events_In_List": 0,
                "Listed_Events": full_events_list,
                "Report_Event_Recorder_Events": full_events_list
            })
        else:
            self.cb.update_params({
                "Report_Event_Recorder_Ignore_Events_In_List": 0,
                "Listed_Events": full_events_list,
                "Report_Event_Recorder_Events": [],
                "Report_Event_Recorder": 1
            })


        # Spatial reporting:
        self.cb.update_params({
            'Enable_Spatial_Output': 1,  # turn on spatial reporting
            # 'Spatial_Output_Channels': ['Infectious_Vectors', 'Adult_Vectors', 'New_Infections', 'Population',
            #                             'Prevalence',
            #                             'New_Diagnostic_Prevalence', 'Daily_EIR', 'New_Clinical_Cases',
            #                             'Human_Infectious_Reservoir', 'Daily_Bites_Per_Human',
            #                             'Land_Temperature',
            #                             'Relative_Humidity', 'Rainfall', 'Air_Temperature']
            'Spatial_Output_Channels': ['Adult_Vectors', 'Population','Prevalence',
                                        'True_Prevalence', 'Daily_EIR', 'Daily_Bites_Per_Human']
        })

        add_summary_report(self.cb)

        # Basic health-seeking

    #################################################################################################
    # REGION-SPECIFIC PARAMETERS

    def africa_setup(self):
        # Use infection/immunity parameters optimized for Africa (courtesy of Prashanth):
        self.cb.update_params({'Antigen_Switch_Rate': pow(10, -9.116590124),
                               'Base_Gametocyte_Production_Rate': 0.06150582,
                               'Base_Gametocyte_Mosquito_Survival_Rate': 0.002011099,

                               'Falciparum_MSP_Variants': 32,
                               'Falciparum_Nonspecific_Types': 76,
                               'Falciparum_PfEMP1_Variants': 1070,
                               'Gametocyte_Stage_Survival_Rate': 0.588569307,

                               'MSP1_Merozoite_Kill_Fraction': 0.511735322,
                               'Max_Individual_Infections': 3,
                               'Nonspecific_Antigenicity_Factor': 0.415111634,
                               })

    #################################################################################################
    # INTERVENTIONS (ADDITIONS TO CAMPAIGN FILE)

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
                                     'seek': 1, 'rate': 0.5},
                                    {'trigger': 'NewSevereCase',
                                     'coverage': 0.1,
                                     'agemin': 5,
                                     'agemax': 100,
                                     'seek': 1,
                                     'rate': 0.5}],
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
        healthseek_events = healthseek_events[np.in1d(healthseek_events['grid_cell'],self.desired_cells)]
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
                                         'seek': 1, 'rate': 0.5},
                                        {'trigger': 'NewSevereCase',
                                         'coverage': float(healthseek_events['cov_severe_adult'][hs]), 'agemin': 5,
                                         'agemax': 100,
                                         'seek': 1,
                                         'rate': 0.5}],
                               drug=['Artemether', 'Lumefantrine'],
                               dosing='FullTreatmentNewDetectionTech',
                               nodes={"class": "NodeSetNodeList", "Node_List": node_list},
                               duration=float(healthseek_events['duration'][hs]))

    def implement_interventions(self, cb, include_itn, include_irs, include_msat, include_mda, include_stepd):
        start_date = "{}-01-01".format(self.start_year)  # Day 1 of simulation
        date_format = "%Y-%m-%d"
        sim_duration = 365 * self.sim_length_years  # length in days
        self.cb.params['Simulation_Duration'] = sim_duration

        # Prevent DTK from spitting out too many messages
        self.cb.params['logLevel_JsonConfigurable'] = "WARNING"
        self.cb.params['Disable_IP_Whitelist'] = 1

        # Event information files
        if not self.itn_fn:
            if include_itn:
                print("WARNING: CANNOT IMPLEMENT ITN INTERVENTION WITHOUT itn_fn SUPPLIED!")
            include_itn = False
        else:
            itn_events = pd.read_csv(self.itn_fn)
            itn_events['simday'] = [convert_to_day(x, start_date, date_format) for x in itn_events.fulldate]
            itn_events = itn_events[np.in1d(itn_events['grid_cell'], self.desired_cells)]
            itn_events.reset_index(inplace=True)

        if not self.irs_fn:
            if include_irs:
                print("WARNING: CANNOT IMPLEMENT IRS INTERVENTION WITHOUT irs_fn SUPPLIED!")
            include_irs = False
        else:
            irs_events = pd.read_csv(self.irs_fn)
            irs_events['simday'] = [convert_to_day(x, start_date, date_format) for x in irs_events.fulldate]
            irs_events = irs_events[np.in1d(irs_events['grid_cell'], self.desired_cells)]
            irs_events.reset_index(inplace=True)

        if not self.msat_fn:
            if include_msat:
                print("WARNING: CANNOT IMPLEMENT MSAT INTERVENTION WITHOUT msat_fn SUPPLIED!")
            include_msat = False
        else:
            msat_events = pd.read_csv(self.msat_fn)
            msat_events['simday'] = [convert_to_day(x, start_date, date_format) for x in msat_events.fulldate]
            msat_events = msat_events[np.in1d(msat_events['grid_cell'], self.desired_cells)]
            msat_events.reset_index(inplace=True)

        if not self.mda_fn:
            if include_mda:
                print("WARNING: CANNOT IMPLEMENT MDA INTERVENTION WITHOUT mda_fn SUPPLIED!")
            include_mda = False
        else:
            mda_events = pd.read_csv(self.mda_fn)
            mda_events['simday'] = [convert_to_day(x, start_date, date_format) for x in mda_events.fulldate]
            mda_events = mda_events[np.in1d(mda_events['grid_cell'], self.desired_cells)]
            mda_events.reset_index(inplace=True)

        if not self.stepd_fn:
            if include_stepd:
                print("WARNING: CANNOT IMPLEMENT STEPD INTERVENTION WITHOUT stepd_fn SUPPLIED!")
            include_stepd = False
        else:
            stepd_events = pd.read_csv(self.stepd_fn)
            stepd_events['simday'] = [convert_to_day(x, start_date, date_format) for x in stepd_events.fulldate]
            stepd_events = stepd_events[np.in1d(stepd_events['grid_cell'], self.desired_cells)]
            stepd_events.reset_index(inplace=True)


        # Get grid cell to node ID lookup table:
        nodeid_lookup,pop_lookup = generate_lookup(self.demo_fp_full)


        if include_itn:
            for itn in range(len(itn_events)):

                # Add non-birth nets
                add_ITN_age_season(cb, start=float(itn_events['simday'][itn]),
                                   age_dep={'youth_cov': float(itn_events['age_cov'][itn]), 'youth_min_age': 5,
                                            'youth_max_age': 20},
                                   coverage_all=float(itn_events['cov_all'][itn]),
                                   as_birth=False,
                                   seasonal_dep={'min_cov': float(itn_events['min_season_cov'][itn]), 'max_day': 1},
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
                                      # dosing='SingleDose',
                                      nodes=[nodeid_lookup[msat_events['grid_cell'][msat]]])

        if include_mda:
            for mda in range(len(mda_events)):
                    add_drug_campaign(cb, campaign_type='MDA', drug_code='DP',
                                      start_days=[float(mda_events['simday'][mda])],
                                      coverage=float(mda_events['cov_all'][mda]), repetitions=1, interval=60,
                                      # dosing='SingleDose',
                                      nodes=[nodeid_lookup[mda_events['grid_cell'][mda]]])

        if include_stepd:
            rcd_people_num = 10 #FIXME should not be hardcoded?
            for sd in range(len(stepd_events)):
                cov = np.min([1.,float(rcd_people_num) / float(pop_lookup[stepd_events['grid_cell'][sd]])])
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


    #################################################################################################
    # ONCE CB IS BUILT, FUNCTIONS FOR WHAT TO DO WITH IT

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

    def submit_experiment(self,num_seeds=1,
                          intervention_sweep=False,
                          migration_sweep=False,
                          vector_migration_sweep=False,
                          simple_intervention_sweep=False,
                          custom_name=None):

        # Implement the actual (not dummy) baseline healthseeking
        self.implement_baseline_healthseeking()


        modlists = []

        if num_seeds > 1:
            new_modlist = [ModFn(DTKConfigBuilder.set_param, 'Run_Number', seed) for seed in range(num_seeds)]
            modlists.append(new_modlist)

        if migration_sweep:
            new_modlist = [ModFn(DTKConfigBuilder.set_param, 'x_Local_Migration', x) for x in [0.5,1,5,10]]
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

        builder = ModBuilder.from_combos(*modlists)

        run_name = self.exp_name
        if custom_name:
            run_name = custom_name


        # SetupParser.init()
        # SetupParser.set("HPC","priority","Normal")
        exp_manager = ExperimentManagerFactory.init()
        exp_manager.run_simulations(config_builder=self.cb, exp_name=run_name, exp_builder=builder)
        return self.cb

    def return_cb_with_interventions(self,
                                     include_itn=False,
                                     include_irs=False,
                                     include_msat=False,
                                     include_mda=False,
                                     include_stepd=False):
        self.implement_baseline_healthseeking()
        self.implement_interventions(self.cb, include_itn, include_irs, include_msat, include_mda, include_stepd)
        return self.cb



class GriddedInputFilesCreator:
    def __init__(self,
                 base,
                 exp_name,
                 desired_cells,
                 cb,
                 grid_pop_csv_file,
                 region="Zambia",
                 start_year=2001,
                 sim_length_years=19,
                 immunity_mode="naive",
                 larval_param_func=None):

        self.base = base
        self.exp_name = exp_name
        self.exp_base = base + 'data/COMPS_experiments/{}/'.format(exp_name)
        self.demo_fp_full = os.path.join(self.exp_base, "Demographics/demo.json")

        self.desired_cells = desired_cells
        self.cb = cb
        self.grid_pop_csv_file = grid_pop_csv_file
        self.immunity_mode = immunity_mode
        self.region = region
        self.start_year = start_year
        self.sim_length_years = sim_length_years
        self.larval_param_func = larval_param_func

        self.gravity_migr_params = np.array([7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00])

        self.file_setup()



    def file_setup(self,generate_immunity_file=True,generate_demographics_file=True,generate_climate_files=True,generate_migration_files=True):

        if generate_demographics_file:
            print("Generating demographics file...")
            self.gen_demo_file()

        # if self.immunity_mode != "naive" and generate_immunity_file:
        #     print("Generating immunity files...")
        #     if self.immunity_mode == "uniform":
        #         self.get_immunity_file_from_single_node()
        #     elif self.immunity_mode == "milen":
        #         self.gen_immunity_file_from_milen_clusters()

        if generate_migration_files:
            print("Generating migration files...")
            self.gen_migration_files()

        if generate_climate_files:
            print("Generating climate files...")
            SetupParser.init()
            [cg_start_year,cg_duration] = safe_start_year_duration_for_climate_generator(self.start_year,self.sim_length_years)

            print("Start year given to climate generator = {}".format(cg_start_year))
            print("Duration given to climate generator = {}".format(cg_duration))

            cg = ClimateGenerator(self.demo_fp_full,
                                  self.exp_base + 'Logs/climate_wo.json',
                                  self.exp_base + 'Climate/',
                                  start_year = str(cg_start_year),
                                  num_years = str(cg_duration),
                                  climate_project = "IDM-{}".format(self.region),
                                  resolution=str(0))

            cg.generate_climate_files()


    def gen_demo_file(self):
        dg = CatchmentDemographicsGenerator.from_file_subset(self.cb, self.grid_pop_csv_file, self.desired_cells)
        demo_dict = dg.generate_demographics()

        # Add catchment name to demographics file metadata:
        # demo_dict["Metadata"]["Catchment"] = self.catch

        # Add larval habitat parameters to demographics file:
        # if self.larval_params_mode != "calibrate":
        demo_dict = self.add_larval_habitats_to_demo(demo_dict)

        # if self.larval_params:
        #     temp_h = self.larval_params['temp_h']
        #     linear_h = self.larval_params['linear_h']
        #     # demo_fp = self.exp_base + "Demographics/demo_temp{}_linear{}.json".format(int(temp_h),int(linear_h))
        #     demo_fp = self.exp_base + "Demographics/demo.json"
        # else:
        demo_fp = self.exp_base + "Demographics/demo.json"

        demo_f = open(demo_fp, 'w+')
        json.dump(demo_dict, demo_f, indent=4)
        demo_f.close()

    def add_larval_habitat_multiplier_to_node(self,node_item, larval_param_dict_this_node):
        calib_single_node_pop = 1000.
        pop_multiplier = float(node_item['NodeAttributes']['InitialPopulation']) / (calib_single_node_pop)

        # Copy the larval param dict handed to this node
        node_item['NodeAttributes']['LarvalHabitatMultiplier'] = larval_param_dict_this_node.copy()

        # Then scale each entry in the dictionary by the population multiplier
        for key in node_item['NodeAttributes']['LarvalHabitatMultiplier'].keys():
            node_item['NodeAttributes']['LarvalHabitatMultiplier'][key] *= pop_multiplier


    def larval_params_func_for_calibration(self, grid_cells):
        return {"CONSTANT": np.ones_like(grid_cells),
                "TEMPORARY_RAINFALL": np.ones_like(grid_cells),
                "LINEAR_SPLINE": np.ones_like(grid_cells),
                "WATER_VEGETATION": np.ones_like(grid_cells)}



    def add_larval_habitats_to_demo(self, demo_dict):
        # Add larval habitat multipliers to demographics file
        # Uses self.larval_params_func

        # This function takes as input grid_cells, and returns a dictionary:
        # {"CONSTANT": [1,2,1,5,...],
        #  "WATER_VEGETATION": etc., }
        # where the list is the multiplier for each node.

        if not self.larval_param_func:
            # Default function, if none is supplied, is the one used for calibration
            self.larval_param_func = self.larval_params_func_for_calibration

        larval_param_multiplier_dict = self.larval_param_func(self.desired_cells)

        larval_param_name_list = list(larval_param_multiplier_dict)
        n_params = len(larval_param_name_list)

        ni = 0
        for node_item in demo_dict['Nodes']:
            larval_param_dict_this_node = {}

            for jj in range(n_params):
                lp_name = larval_param_name_list[jj]
                larval_param_dict_this_node[lp_name] = larval_param_multiplier_dict[lp_name][jj]

            self.add_larval_habitat_multiplier_to_node(node_item, larval_param_dict_this_node)

            ni += 1

        return demo_dict


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
                                                 self.exp_base +'Migration/local_migration.bin.json')



class CatchmentDemographicsGenerator(DemographicsGenerator):

    @classmethod
    def from_file_subset(cls, cb, population_input_file, desired_cells,demographics_type='static', res_in_arcsec=30,
                         update_demographics=None, default_pop=1000):

        nodes_list = list()
        with open(population_input_file, 'r') as pop_csv:
            reader = csv.DictReader(pop_csv)
            for row in reader:
                # Latitude
                if not 'lat' in row: raise ValueError('Column lat is required in input population file.')
                lat = float(row['lat'])

                # Longitude
                if not 'lon' in row: raise ValueError('Column lon is required in input population file.')
                lon = float(row['lon'])

                # Node label
                res_in_deg = cls.arcsec_to_deg(res_in_arcsec)
                node_label = row['node_label'] # if 'node_label' in row else nodeid_from_lat_lon(lat, lon, res_in_deg)

                # Population
                pop = int(float(row['pop'])) if 'pop' in row else default_pop

                if int(node_label) in desired_cells:
                    # Append the newly created node to the list
                    nodes_list.append(Node(lat, lon, pop, node_label))

        return cls(cb, nodes_list, demographics_type, res_in_arcsec, update_demographics, default_pop)
