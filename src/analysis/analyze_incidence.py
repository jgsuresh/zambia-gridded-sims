from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from relative_time import *
from simtools.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
import numpy as np
import pandas as pd

from gridded_sim_general import *
from RDT_map_from_data import *
import datetime
from datetime import date


class RDTIncidenceAnalyzer(BaseAnalyzer):

    # filenames = ['output\SpatialReport_Population.bin',
    #              'output\SpatialReport_Prevalence.bin',
    #              'output\SpatialReport_New_Diagnostic_Prevalence.bin',
    #              'Assets\Demographics\demo.json']
    filenames = ['output/SpatialReport_New_Clinical_Cases.bin',
                 'output/ReportEventRecorder.csv',
                 'Assets/Demographics/demo.json']

    def __init__(self):
        super(RDTIncidenceAnalyzer, self).__init__()
        self.my_data = {}
        self.metadata = {}

        self.cases_by_node = {}
        self.treated_by_node = {}
        self.catch = {}
        self.demo_file = {}
        self.node_ids = {}

    def filter(self, sim_metadata):
        return sim_metadata['Run_Number'] == 0

    def apply(self, parser):
        exp_name = parser.experiment.exp_name
        self.catch = exp_name.split('_')[0] # Assumes the experiment name is "CATCHNAME_otherstuff"

        #fixme: In the future, use the parser to access the demographics file from Assets directly.
        base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
        self.demo_file = base + "data/COMPS_experiments/{}_full/Demographics/demo.json".format(self.catch)

        case_data = parser.raw_data[self.filenames[0]]
        recorder = parser.raw_data[self.filenames[1]]
        # demo_data = parser.raw_data[self.filenames[2]]

        self.node_ids = case_data['nodeids']
        self.n_tstep = case_data['n_tstep']
        self.n_nodes = case_data['n_nodes']

        self.cases_by_node = {}
        self.treated_by_node = {}
        # Initialize node arrays:
        for j in range(self.n_nodes):
            self.cases_by_node[j] = np.zeros(self.n_tstep)
            self.treated_by_node[j] = np.zeros(self.n_tstep)

        # Collect per-node data from custom reporter for cases:
        for i in range(self.n_tstep):
            case_timestep_data = case_data['data'][i]
            for j in range(self.n_nodes):
                self.cases_by_node[j][i] = case_timestep_data[j]

        # Collect per-node data from event recorder for number treated:
        self.weekly_treated_sim = aggregate_events_in_recorder(recorder,
                                                               "Received_Treatment",
                                                               aggregate_type="week",
                                                               start_date="2007-01-01")

        self.monthly_treated_sim = aggregate_events_in_recorder(recorder,
                                                                "Received_Treatment",
                                                                aggregate_type="month",
                                                                start_date="2007-01-01")

    def finalize(self):
        print("")

    def plot(self):

        def check_alternate_catch_names(test_str,catch):
            pass

        def load_DHIS_HF_and_CHW_data(base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'): #merge_clinical_and_rdt_for_HF=False):
            print("Loading DHIS data...")
            obs_df = pd.read_csv(base + "data/incidence/Oct2017_AllSouthern.csv",
                                 usecols=['dataelement', 'period', 'orgunit', 'value'])
            keep_data_elements = ['Clinical malaria cases', 'Passive Number Positive', 'RDT positive cases']
            keep_data = obs_df.apply(lambda x: x['dataelement'] in keep_data_elements, axis=1)
            obs_df = obs_df[keep_data]

            in_catch = obs_df.apply(lambda x: (("so {}".format(self.catch.title()) in x['orgunit']) or
                                               ("so {}".format(self.catch) in x['orgunit'])), axis=1)
            # if self.catch == "nyanga chaamwe":
            #     in_catch = obs_df.apply(lambda x: ("so NyangaC".format(self.catch.capitalize()) in x['orgunit'],axis=1)
            catch_df = obs_df[in_catch]
            # fixme Does not take into account the possibility of name variations (e.g. Lukonde vs Lukande).  Could accomplish this with an "acceptable name variations" dictionary
            print("Done loading DHIS data.")

            in_HF = catch_df.apply(lambda x: "W" in x['period'], axis=1)
            HF_df = catch_df[in_HF]
            CHW_df = catch_df[np.logical_not(in_HF)]

            # if merge_clinical_and_rdt_for_HF:
            #     HF_df["Full clinical cases"] = HF_df["Clinical malaria cases"] + HF_df["RDT positive cases"]

            return [HF_df,CHW_df]

        def add_weeks_and_years_to_HF_df(HF_df):
            HF_df['year_num'] = HF_df.apply(lambda x: np.int32(x['period'].split("W")[0]), axis=1)
            HF_df['week_num'] = HF_df.apply(lambda x: np.int32(x['period'].split("W")[1]), axis=1)
            return HF_df

        def add_months_and_years_to_CHW_df(CHW_df):
            CHW_df['year_num'] = CHW_df.apply(lambda x: np.int32(x['period'][:4]), axis=1)
            CHW_df['month_num'] = CHW_df.apply(lambda x: np.int32(x['period'][4:]), axis=1)
            return CHW_df

        def generate_week_index():
            # Generate week index from start to end:
            start_year = np.int32(start_date.split('-')[0])
            sim_duration = self.n_tstep
            sim_duration_years = self.n_tstep/365 # Assumes that the simulation is running for some number of years

            sim_duration_weeks = sim_duration_years * 52
            week_counter = np.arange(sim_duration_weeks)
            week_num = np.tile(np.arange(52)+1,sim_duration_years)
            year_num = np.repeat(np.arange(start_year,start_year + sim_duration_years),52)

            # Make a dataframe of this stuff that we can add to:
            weeks_df = pd.DataFrame({
                "week_counter": week_counter,
                "week_num": week_num,
                "year_num": year_num
            })

            return [week_counter, week_num, year_num, weeks_df]

        def aggregate_HF_cases_sim():
            # Add up all weekly cases across all nodes in the simulation
            test = self.weekly_treated_sim.groupby(['year','week']).sum()['count']
            test = test.reset_index()
            return test

        # def convert_week_to_date(week_num,year_num):
        #     date_string = "{}-W{}".format(year_num,week_num)
        #     return datetime.strptime(date_string + '-0', "%Y-W%W-%w")
        def convert_week_to_date(week_num,year_num):
            date_string = "{}-W{}".format(year_num,week_num)
            date_obj = datetime.datetime.strptime(date_string + '-0', "%Y-W%W-%w")
            date_str = "{}-{}-{}".format(date_obj.year,str(date_obj.month).zfill(2),str(date_obj.day).zfill(2))
            return date_str

        def convert_month_to_date(month_num,year_num):
            date_str = "{}-{}-01".format(year_num,str(month_num).zfill(2))
            return date_str

        # ============================================================================================================
        # ============================================================================================================

        import matplotlib.pyplot as plt
        from matplotlib import cm
        import matplotlib.dates as mdates
        import seaborn as sns
        sns.set_style("darkgrid")

        start_date = "2007-01-01"  # Day 1 of simulation
        date_format = "%Y-%m-%d"
        base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'

        plot_HF = False
        plot_CHW = True

        plot_sim = True
        plot_data = True

        if plot_data:
            [HF_df, CHW_df] = load_DHIS_HF_and_CHW_data(base=base)  # merge_clinical_and_rdt_for_HF=True)
            HF_df = add_weeks_and_years_to_HF_df(HF_df)
            CHW_df = add_months_and_years_to_CHW_df(CHW_df)

        # HEALTH FACILITIES:
        if plot_HF:
            # Look up DHIS incidence data, which will be used for comparison:

            [week_counter, week_num, year_num, weeks_df] = generate_week_index()
            # print weeks_df

            if plot_data:
                # Merge in the observational data to the combined dataframe
                weeks_df = weeks_df.merge(HF_df[HF_df['dataelement']=='RDT positive cases'], how='left', left_on=['year_num','week_num'], right_on=['year_num','week_num'])
                # weeks_df = weeks_df.merge(HF_df[HF_df['dataelement'] == 'Full clinical cases'], how='left', left_on=['year_num', 'week_num'], right_on=['year_num', 'week_num'])
                weeks_df = weeks_df.drop(["dataelement","period"],axis=1)
                weeks_df = weeks_df.rename(columns={"value": "Data_RDT"})
                # weeks_df = weeks_df.rename(columns={"value": "Data"})

                # For Health Facilities, desired case counts is taken to be sum of "Clinical malaria cases" and "RDT positive cases"
                weeks_df = weeks_df.merge(HF_df[HF_df['dataelement'] == 'Clinical malaria cases'], how='left',left_on=['year_num', 'week_num'], right_on=['year_num', 'week_num'])
                weeks_df = weeks_df.rename(columns={"value": "Data_Clinical"})
                weeks_df["Data"] = weeks_df["Data_RDT"] + weeks_df["Data_Clinical"]
                weeks_df = weeks_df.drop(["dataelement", "period"], axis=1)


            # Merge in the simulation data to this dataframe:
            self.weekly_treated_HF_sim = aggregate_HF_cases_sim()
            weeks_df = weeks_df.merge(self.weekly_treated_HF_sim, how='left',left_on=['year_num','week_num'], right_on=['year','week'])
            weeks_df = weeks_df.rename(columns={"count": "Sim"})

            # Missing data from observations should be shown as missing.  But missing data from sims are zeros.
            weeks_df['Sim'] = weeks_df['Sim'].fillna(0)


            # For plotting purposes: get weeks in terms of actual dates:
            foo = mdates.strpdate2num(date_format)
            weekdates_mdates = np.array([])
            for i in week_counter:
                hold = convert_week_to_date(week_num[i],year_num[i])
                weekdates_mdates = np.append(weekdates_mdates,foo(hold))

            # Plot!!!
            plt.figure(figsize=(15,5))
            if plot_sim:
                # plt.plot(weeks_df['week_counter'],weeks_df['Sim'],label="Sim Incidence")
                plt.plot_date(weekdates_mdates, weeks_df['Sim'], label="Sim Incidence",
                              linestyle='-',marker=',',lw=2.5)
            if plot_data:
                # plt.plot(weeks_df['week_counter'],weeks_df['Data'],label="Data Incidence")
                plt.plot_date(weekdates_mdates, weeks_df['Data'], label="Data Incidence (sum of RDT+ and Clinical)",
                              linestyle='-',marker=',',lw=2.5)
                plt.plot_date(weekdates_mdates, weeks_df['Data_RDT'], label="DHIS: RDT positive cases]",
                              linestyle='dotted', marker=',')
                plt.plot_date(weekdates_mdates, weeks_df['Data_Clinical'], label="DHIS: Clinical malaria cases]",
                              linestyle='dotted', marker=',')

            plt.xlabel("Week Number")
            plt.ylabel("HF weekly cases")
            plt.legend()
            plt.title(self.catch.capitalize())
            plt.show()
            # plt.savefig(base + "data/figs/incidence/{}_HF.png".format(self.catch))


        if plot_CHW:
            def get_list_of_CHWs(base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
                # Open Caitlin's file
                lookup_df = pd.read_csv(base + "data/incidence/CHW_grid_lookup.csv")
                in_catch = lookup_df.apply(lambda x: "so {}".format(self.catch.title()) in x['closest.chw.name'],axis=1)
                catch_df = lookup_df[in_catch]
                return list(set(catch_df['closest.chw.name']))

            def get_grid_cells_for_CHW(chw_name,base='C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'):
                lookup_df = pd.read_csv(base + "data/incidence/CHW_grid_lookup.csv")

                this_CHW = lookup_df["closest.chw.name"] == chw_name
                # Sometimes, the CHW has "chw_name" with A/B/C etc added afterwards.  These variations need to be added together.
                #fixme Actually this only needs to be taken into account on the data side, not the simulation side

                # this_CHW = lookup_df.apply(lambda x: chw_name in x["closest.chw.name"],axis=1)

                # I believe that this should work even if there are more than one CHW that we are putting together:
                # if len(list(set(lookup_df[this_CHW]["closest.chw.name"]))) == 1:
                return lookup_df['loc.id'][this_CHW]

            def get_node_ids_for_CHW(chw_name):
                grid_cells = get_grid_cells_for_CHW(chw_name)
                dtk_node_ids = convert_from_grid_cells_to_dtk_node_ids_using_demo(grid_cells,self.demo_file)
                return dtk_node_ids

            def generate_month_index():
                # Generate week index from start to end:
                start_year = np.int32(start_date.split('-')[0])
                sim_duration = self.n_tstep
                sim_duration_years = self.n_tstep / 365  # Assumes that the simulation is running for some number of years

                sim_duration_months = sim_duration_years * 12
                month_counter = np.arange(sim_duration_months)
                month_num = np.tile(np.arange(12) + 1, sim_duration_years)
                year_num = np.repeat(np.arange(start_year, start_year + sim_duration_years), 12)

                # Make a dataframe of this stuff that we can add to:
                months_df = pd.DataFrame({
                    "month_counter": month_counter,
                    "month_num": month_num,
                    "year_num": year_num
                })

                return [month_counter, month_num, year_num, months_df]

            def aggregate_CHW_cases_sim(chw_name):
                # Add up all montly cases across CHW-specific nodes in the simulation
                chw_nodeset = get_node_ids_for_CHW(chw_name)
                this_chw = self.monthly_treated_sim.apply(lambda x: x['Node_ID'] in chw_nodeset,axis=1)
                test = self.monthly_treated_sim[this_chw].groupby(['year', 'month']).sum()['count']
                test = test.reset_index()
                return test

            def aggregate_CHW_df_data(CHW_df,chw_name,dataelement='Passive Number Positive'):
                # Get compiled dataframe for CHW's in data.  This can sometimes be complicated by people like "Bob A" and "Bob B"
                CHW_df = CHW_df[CHW_df['dataelement']==dataelement]

                this_CHW = CHW_df.apply(lambda x: chw_name in x["orgunit"], axis=1)
                if np.sum(this_CHW) == 0:
                    print("ZERO instances found of ",chw_name)
                this_CHW_df = CHW_df[this_CHW]

                # Check how many CHWs there are that fit this name:
                if len(list(set(this_CHW_df['orgunit']))) == 1:
                    this_CHW_df = this_CHW_df.drop(["dataelement","period"],axis=1)
                    return this_CHW_df
                else:
                    # There are more than 1.  Will need to aggregate:
                    # Note that this process automatically drops the other unnecessary columns
                    agg_CHW_df = this_CHW_df.groupby(['year_num','month_num']).sum()['value']
                    agg_CHW_df = agg_CHW_df.reset_index()
                    return agg_CHW_df

            # Make dataframe structure to house data from both observations and simulations:
            [month_counter, month_num, year_num, months_df] = generate_month_index()

            chw_list = get_list_of_CHWs()
            print("chw_list: ",chw_list)
            num_chw = len(chw_list)
            plt.figure(figsize=(10,10))

            chw_counter = 1
            for chw_name in chw_list:
                # chw_name = "so Chabbobboma RHC Syambabala"
                print("Working on CHW {}...".format(chw_name))

                # Merge in observations:
                # print CHW_df
                # this_CHW_df = CHW_df[CHW_df['orgunit'] == chw_name]
                this_CHW_df = aggregate_CHW_df_data(CHW_df,chw_name)
                months_this_chw_df = months_df.merge(this_CHW_df, how='left', left_on=['year_num','month_num'], right_on=['year_num','month_num'])
                months_this_chw_df = months_this_chw_df.rename(columns={"value": "Data"})

                # Merge in simulation results:
                monthly_treated_this_CHW_sim = aggregate_CHW_cases_sim(chw_name)
                months_this_chw_df = months_this_chw_df.merge(monthly_treated_this_CHW_sim, how='left',left_on=['year_num','month_num'], right_on=['year','month'])
                months_this_chw_df = months_this_chw_df.rename(columns={"count": "Sim"})

                # Missing data from observations should be shown as missing.  But missing data from sims are zeros.
                months_this_chw_df['Sim'] = months_this_chw_df['Sim'].fillna(0)

                # For plotting purposes: get months in terms of actual dates:
                foo = mdates.strpdate2num(date_format)
                monthdates_mdates = np.array([])
                for i in month_counter:
                    hold = convert_month_to_date(month_num[i],year_num[i])
                    monthdates_mdates = np.append(monthdates_mdates,foo(hold))

                # Plot:
                ax = plt.subplot(num_chw, 1, chw_counter)
                # ax.plot(month_counter, months_this_chw_df['Sim'], label="Sim Incidence", linestyle='-')
                # ax.plot(month_counter, months_this_chw_df['Data'], label="Data Incidence", linestyle='-')
                ax.plot_date(monthdates_mdates, months_this_chw_df['Sim'], label="Sim Incidence", linestyle='-',marker=',')
                ax.plot_date(monthdates_mdates, months_this_chw_df['Data'], label="Data Incidence", linestyle='-',marker=',')
                ax.legend()
                ax.set_xlabel("Months")
                ax.set_ylabel("Monthly case counts")
                ax.set_title(chw_name)

                chw_counter += 1

            plt.tight_layout()
            # plt.show()
            plt.savefig(base + "data/figs/incidence/{}_CHW.png".format(self.catch))

if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager.AnalyzeManager()

    # am.add_experiment(retrieve_experiment("43cac760-cbd6-e711-9414-f0921c16b9e5")) # bbondo
    # am.add_experiment(retrieve_experiment("a31b516a-cbd6-e711-9414-f0921c16b9e5"))  # chabbobboma
    # am.add_experiment(retrieve_experiment("1ecdf372-cbd6-e711-9414-f0921c16b9e5")) # chisanga
    # am.add_experiment(retrieve_experiment("957e6159-32d6-e711-9414-f0921c16b9e5")) # chiyabi
    # am.add_experiment(retrieve_experiment("9669907b-cbd6-e711-9414-f0921c16b9e5"))  # luumbo
    # am.add_experiment(retrieve_experiment("fbe40809-ccd6-e711-9414-f0921c16b9e5"))  # munyumbwe
    # am.add_experiment(retrieve_experiment("8aadd6a0-cbd6-e711-9414-f0921c16b9e5"))  # nyanga chaamwe
    # am.add_experiment(retrieve_experiment("d18a9aa8-cbd6-e711-9414-f0921c16b9e5"))  # sinafala
    am.add_experiment(retrieve_experiment("d28a9aa8-cbd6-e711-9414-f0921c16b9e5"))  # sinamalima


    am.add_analyzer(RDTIncidenceAnalyzer())
    am.analyze()