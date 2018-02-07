import json

import pandas as pd
import numpy as np

from calibtool.analyzers.BaseCalibrationAnalyzer import BaseCalibrationAnalyzer

from calibtool.analyzers.BaseComparisonAnalyzer import BaseComparisonAnalyzer

from gridded_sim_general import *
from relative_time import *


class CalibMozambiqueAnalyzer(BaseCalibrationAnalyzer):
    filenames = ['output/SpatialReport_Population.bin', 'output/SpatialReport_New_Diagnostic_Prevalence.bin',
                 'Assets/Demographics/demo.json']
    # filenames = ['output/InsetChart.json']

    def __init__(self, site):
        super().__init__(site)
        self.start_date = "1994-01-01"
        self.base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
        self.error = {}


    # def apply(self, parser):
    #     sample_index = parser.sim_data.get('__sample_index__')
    #
    #     df = pd.DataFrame()
    #     df["RDT"] = np.array([1.,2.,3.,4.,5.,6.])
    #     df.sample = sample_index
    #     df.sim_id = parser.sim_id
    #     return df

    # def filter(self, sim_metadata):
    #     return True

    def apply(self, parser):
        sample_index = parser.sim_data.get('__sample_index__')
        run_number = parser.sim_data['Run_Number']

        pop_data = parser.raw_data[self.filenames[0]]
        RDT_prev_data = parser.raw_data[self.filenames[1]]
        demo = parser.raw_data[self.filenames[2]]

        node_ids = pop_data['nodeids']
        n_tstep = pop_data['n_tstep']
        n_nodes = pop_data['n_nodes']

        # Get initial population of nodes:
        pop_init = np.zeros(n_nodes)
        for ni in range(n_nodes):
            pop_init[ni] = pop_data['data'][0][ni]

        # Get list of grid cells that correspond to these nodes
        grid_cells = convert_from_dtk_node_ids_to_grid_cells_using_demo(node_ids,demo)
        # From list of grid cell IDs, get round dates, populations (which will be used for weight), and prevalence
        grid_prev_df = get_RDT_ref_data_for_grid_cells(grid_cells, path_from_base="data/mozambique/grid_prevalence_with_dates.csv")

        # Get corresponding prevalence values from simulation (match by grid cell ID, same date [minus 1 day])
        node_grid_df = pd.DataFrame({
            "node": node_ids,
            "node_index": np.arange(len(node_ids)),
            "grid_cell": grid_cells
        })
        prev_comparison_df = grid_prev_df.merge(node_grid_df,how='left',left_on='grid_cell',right_on='grid_cell')

        n_obs = len(prev_comparison_df)

        day_num = np.zeros(n_obs)

        # Convert observed round dates to simulation day numbers:
        date_list = list(prev_comparison_df["date"])
        ii = 0
        for date in date_list:
            day_num[ii] = convert_to_day_365(date,self.start_date)
            ii += 1

        # Account for instantaneous MDA drops by shifting the day numbers back by 1:
        day_num = day_num-1

        # Only do comparison when we actually have both sim and observations:
        comparison_range = np.logical_and(day_num >= 0, day_num < n_tstep)
        # print("Total days in comparison: ", np.sum(comparison_range))

        if True: # this is correct way to do it, for future
            day_num = day_num[comparison_range]
        if False: # for debugging purposes only
            comparison_range[:9] = np.ones(9, dtype=bool)
            day_num = np.arange(9)

        prev_comparison_df = prev_comparison_df[comparison_range]
        prev_comparison_df = prev_comparison_df.reset_index()

        node_index = np.array(prev_comparison_df["node_index"])
        n_comparisons = np.sum(comparison_range)
        prev_sim = np.zeros(n_comparisons)

        # Get corresponding prevalence values from simulation
        for ii in range(n_comparisons):
            dn = int(day_num[ii])
            ni = node_index[ii]
            # try:
            prev_sim[ii] = RDT_prev_data['data'][dn][ni]
            # except:
            #     print("FAILED: ",dn, ni)

        # Add the simulation prevalence values to the comparison dataframe
        # prev_comparison_df["prev_sim"] = pd.Series(prev_sim,index=prev_comparison_df.index)
        # col_name = "prev_sim_sample{}_run{}".format(sample_index,run_number)
        col_name = "prev_sim_sample{}".format(sample_index)
        prev_comparison_df[col_name] = pd.Series(prev_sim, index=prev_comparison_df.index)

        prev_comparison_df.sample = sample_index
        prev_comparison_df.sim_id = parser.sim_id

        parser.prev_comparison_df = prev_comparison_df.copy()
        self.holdme = prev_comparison_df.copy()
        return prev_comparison_df

        # # Compute a pop-weighted MSE here (see green whiteboard)
        # self.error[sample_index] = self._MSE_calc(prev_comparison_df[comparison_range])



        # sim_data_df = pd.DataFrame({
        #     "RDT_sim": prev_sim[comparison_range],
        #     "day_num": day_num[comparison_range],
        #     "pop_init": pop_init[comparison_range] #fixme Not sure if this is correct
        #     "RDT_
        # })
        # sim_data_df.sample = sample_index
        # sim_data_df.sim_id = parser.sim_id
        # return sim_data_df


    def combine(self, parsers):
        """
        Combine the simulation data into a single table for all analyzed simulations.
        """
        # selected = [p.selected_data[id(self)] for p in parsers.values() if id(self) in p.selected_data]

        self.n_samples = len(parsers)
        # Merge all of the dataframes together:

        for sim_id, parser in parsers.items():
            if 'data' in locals():
                data = data.merge(parser.prev_comparison_df.copy(),
                                  how='left',
                                  left_on=['grid_cell', 'date', 'N', 'prev'],
                                  right_on=['grid_cell', 'date', 'N', 'prev'])
            else:
                data = parser.prev_comparison_df.copy()

        sample_col_list = ["prev_sim_sample{}".format(sample_ind) for sample_ind in range(self.n_samples)]
        other_col_list = ["grid_cell","N","date","prev"]
        full_col_list = sample_col_list + other_col_list

        self.data = data[full_col_list].copy()

        # Make sure to drop any NA's that may be lurking in the reference data
        self.data = self.data.dropna()

    def cache(self):
        return self.data

    def plot(self):
        pass

    def plot_comparison(cls, fig, data, **kwargs):
        pass

    def compare(self, sample):
        """
        Assess the result per sample, in this case the likelihood
        comparison between simulation and reference data.
        """
        return -1

    def finalize(self):
        """
        Calculate the output result for each sample.
        """

        result_arr = np.zeros(self.n_samples)
        for i in range(self.n_samples):
            # sample_col_name = "prev_sim_sample{}".format(i)

            # Minus sign very important, because OptimTool is trying to maximize the result
            result_arr[i] = -1 * self.calc_mse_sqrt(self.data,i)

        self.result = pd.Series(result_arr)
        self.result.index.name = "sample"
        print(self.result)


    def calc_mse_sqrt(self, data_df, sample_index):
        sim_data = np.array(data_df["prev_sim_sample{}".format(sample_index)])
        ref_data = np.array(data_df["prev"])
        pops = np.array(data_df["N"])

        mse = np.sum(pops * (sim_data - ref_data)**2.)/np.sum(pops)
        mse_sqrt = np.sqrt(mse)

        return mse_sqrt

    # @classmethod
    # def plot_comparison(cls, fig, data, **kwargs):
    #     ax = fig.gca()
    #
    #     df = pd.DataFrame.from_dict(data, orient='columns')
    #     for iax, (species, group_df) in enumerate(df.groupby('Channel')):
    #         counts = list(group_df['Counts'])
    #         time = list(group_df['Month'])
    #         months = [calendar.month_abbr[i] for i in time]
    #         if kwargs.pop('reference', False):
    #             logger.debug('months: %s', time)
    #             logger.debug('counts: %s', counts)
    #             ax.plot(months, counts, **kwargs)
    #         else:
    #             fmt_str = kwargs.pop('fmt', None)
    #             args = (fmt_str,) if fmt_str else ()
    #             ax.plot(months, counts, *args, **kwargs)
    #     ax.set(xlabel='Months', ylabel=species)



    # def finalize(self):
    #     # self.result = self.data.groupby(level='sample', axis=1).apply(self.compare)
    #     self.result = pd.Series(np.arange(5))
    #
    # def cache(self):
    #     ref = np.array([3.,4.,5.,6.,7.,8.])
    #     data = np.array([1.,2.,3.,4.,5.,6.])
    #
    #     return {"ref": ref, "data": data}
    #
    #
    # def compare(self, sample):
    #     return 5