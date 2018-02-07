import json

import pandas as pd
import numpy as np

from calibtool.analyzers.BaseCalibrationAnalyzer import BaseCalibrationAnalyzer

from calibtool.analyzers.BaseComparisonAnalyzer import BaseComparisonAnalyzer

from gridded_sim_general import *
from relative_time import *


class GriddedRDTLikelihoodAnalyzer(BaseCalibrationAnalyzer):
    filenames = ['output/SpatialReport_Population.bin', 'output/SpatialReport_New_Diagnostic_Prevalence.bin',
                 'Assets/Demographics/demo.json']

    def __init__(self, site):
        super().__init__(site)
        self.start_date = "2010-01-01"
        self.base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
        self.error = {}


    def apply(self, parser):
        """
        Extract data from output data and accumulate in same bins as reference.
        """

        # Load data from simulation
        data = parser.raw_data[self.filenames[0]]

        data = data[365:]
        data['Day'] = data['Time'].apply(lambda x: (x + 1) % 365)
        data = data[['Day', 'Species', 'Population', 'VectorPopulation']]
        data['Vector_per_Human'] = data['VectorPopulation'] / data['Population']
        data = data.groupby(['Day', 'Species'])['Vector_per_Human'].apply(np.mean).reset_index()

        dateparser = lambda x: datetime.datetime.strptime(x, '%j').month
        data['Month'] = data['Day'].apply(lambda x: dateparser(str(x + 1)))
        data = data.groupby(['Month', 'Species'])['Vector_per_Human'].apply(np.mean).reset_index()

        data = data.rename(columns={'Vector_per_Human': 'Counts', 'Species': 'Channel'})
        data = data.sort_values(['Channel', 'Month'])
        data = data.set_index(['Channel', 'Month'])
        channel_data_dict = {}

        for channel in self.site.metadata['species']:

            with thread_lock:  # TODO: re-code following block to ensure thread safety (Issue #758)?

                # Reset multi-index and perform transformations on index columns
                df = data.copy().reset_index()
                df = df.rename(columns={'Counts': channel})
                del df['Channel']

                # Re-bin according to reference and return single-channel Series
                rebinned = aggregate_on_index(df, self.reference.loc(axis=1)[channel].index, keep=[channel])
                channel_data_dict[channel] = rebinned[channel].rename('Counts')

        sim_data = pd.concat(channel_data_dict.values(), keys=channel_data_dict.keys(), names=['Channel'])
        sim_data = pd.DataFrame(sim_data)  # single-column DataFrame for standardized combine/compare pattern
        sim_data.sample = parser.sim_data.get('__sample_index__')
        sim_data.sim_id = parser.sim_id

        return sim_data
    # Prashanth
    """
    def apply(self, parser):
        # Extract data from output data and accumulate in same bins as reference.
        

        # Load data from simulation
        data = parser.raw_data[self.filenames[0]]

        data = data[365:]
        data['Day'] = data['Time'].apply(lambda x: (x + 1) % 365)
        data = data[['Day', 'Species', 'Population', 'VectorPopulation']]
        data['Vector_per_Human'] = data['VectorPopulation'] / data['Population']
        data = data.groupby(['Day', 'Species'])['Vector_per_Human'].apply(np.mean).reset_index()

        dateparser = lambda x: datetime.datetime.strptime(x, '%j').month
        data['Month'] = data['Day'].apply(lambda x: dateparser(str(x + 1)))
        data = data.groupby(['Month', 'Species'])['Vector_per_Human'].apply(np.mean).reset_index()

        data = data.rename(columns={'Vector_per_Human': 'Counts', 'Species': 'Channel'})
        data = data.sort_values(['Channel', 'Month'])
        data = data.set_index(['Channel', 'Month'])
        channel_data_dict = {}

        for channel in self.site.metadata['species']:

            with thread_lock:  # TODO: re-code following block to ensure thread safety (Issue #758)?

                # Reset multi-index and perform transformations on index columns
                df = data.copy().reset_index()
                df = df.rename(columns={'Counts': channel})
                del df['Channel']

                # Re-bin according to reference and return single-channel Series
                rebinned = aggregate_on_index(df, self.reference.loc(axis=1)[channel].index, keep=[channel])
                channel_data_dict[channel] = rebinned[channel].rename('Counts')

        sim_data = pd.concat(channel_data_dict.values(), keys=channel_data_dict.keys(), names=['Channel'])
        sim_data = pd.DataFrame(sim_data)  # single-column DataFrame for standardized combine/compare pattern
        sim_data.sample = parser.sim_data.get('__sample_index__')
        sim_data.sim_id = parser.sim_id

        return sim_data
    """

    def apply(self, parser):
        sample_index = parser.sim_data.get('__sample_index__')

        pop_data = parser.raw_data[self.filenames[0]]
        RDT_prev_data = parser.raw_data[self.filenames[1]]
        demo = parser.raw_data[self.filenames[2]]

        # Parse demographics file for catchment name
        # catch = demo["Metadata"]["Catchment"]
        #fixme Can also parse demo file for initial populations in the future

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
        grid_prev_df = get_RDT_ref_data_for_grid_cells(grid_cells)
        n_obs = len(grid_prev_df)

        # Get corresponding prevalence values from simulation (match by grid cell ID, same date [minus 1 day])
        node_grid_df = pd.DataFrame({
            "node": node_ids,
            "node_index": np.arange(len(node_ids)),
            "grid": grid_cells
        })
        prev_comparison_df = grid_prev_df.merge(node_grid_df,how='left',left_on='grid_cell',right_on='grid')
        node_index = np.array(prev_comparison_df["node_index"])

        prev_sim = np.zeros(len(grid_cells))
        day_num = np.zeros(len(grid_cells))

        # Convert observed round dates to simulation day numbers:
        date_list = list(prev_comparison_df["date"])
        ii = 0
        for date in date_list:
            day_num[ii] = convert_to_date_365(date,self.start_date)
            ii += 1

        # Account for instantaneous MDA drops:
        day_num = day_num-1

        # Get corresponding prevalence values from simulation
        for ii in range(n_obs):
            prev_sim[ii] = RDT_prev_data['data'][day_num[ii]][node_index[ii]]

        # Add the simulation prevalence values to the comparison dataframe
        prev_comparison_df["prev_sim"] = pd.Series(prev_sim,index=prev_comparison_df.index)

        # Only do comparison when we actually have both sim and obsevations:
        comparison_range = np.logical_and(day_num >= 0, day_num < n_tstep)

        # Compute a pop-weighted MSE here (see green whiteboard)
        self.error[sample_index] = self._MSE_calc(prev_comparison_df[comparison_range])



    def _MSE_calc(self,prev_comparison_df):
        w = np.array(prev_comparison_df['pop'])
        diffsq = (np.array(prev_comparison_df['prev'])-np.array(prev_comparison_df['prev_sim']))**2
        weighted_diffsq = np.sum(w*diffsq)/np.sum(w)
        return -np.sqrt(weighted_diffsq) #negative sign because OptimTool tries to maximize

    def combine(self, parsers):
        pass
    # def combine(self, parsers):
    #     # Collect all the data for all the simulations
    #     for p in parsers.values():
    #         self.data.append(p.selected_data[id(self)])
    #
    #     # Sort our data by sample_index
    #     # We need to preserve the order by sample_index
    #     self.data = sorted(self.data, key=lambda k: k['sample_index'])

    def plot(self):
        pass

    def finalize(self):
        sorted_sample_index_list = sorted(self.error) #sorted(self.error) #, key=lambda k: k['sample_index'])

        result_list = []
        for sample_index in sorted_sample_index_list:
            result_list.append(self.error[sample_index])

        # Result needs to be a series where the index is the sample_index and the value the likelihood for this sample
        self.result = pd.Series(result_list, index=sorted_sample_index_list)



    # def cache(self):
    #     return json.dumps(self.data)
