# Create lookup table between grid IDs and node IDs or populations in gridded simulations:

import pandas as pd
import numpy as np
from gen_migr_json import *

def generate_lookup(demo_file_fp):
    df = load_demo(demo_file_fp)

    id_lookup = {}
    pop_lookup = {}
    for i,row in df.iterrows():
        id_lookup[int(row['grid_id'])] = int(row['node_id'])
        pop_lookup[int(row['grid_id'])] = int(row['pop'])
    return [id_lookup,pop_lookup]

def generate_latlong_lookup(demo_file_fp):
    df = load_demo(demo_file_fp)

    lookup = {}
    for i, row in df.iterrows():
        lookup[int(row['node_id'])] = [row['lat'],row['long']]
    return lookup

# Testing:
# if __name__ == "__main__":
#     foo = generate_lookup('C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/experiments/gravity_test_v0/inputs/Demographics/MultiNode/demo.json')
#     print foo