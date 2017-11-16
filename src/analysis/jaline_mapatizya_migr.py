# Compute migration rates for Mapatizya (for Jaline's household model setup)

import pandas as pd
import numpy as np
from geopy.distance import vincenty
import json

def compute_migr_prob(grav_params, ph, pd, d):
    num_trips = grav_params[0] * ph ** grav_params[1] * pd ** grav_params[2] * d ** grav_params[3]
    prob_trip = np.min([1., num_trips / ph])
    return prob_trip


def compute_migr_dict(df, grav_params):
    migr = {}

    p_sum = np.zeros(len(df))
    jj = 0
    for i1, r1 in df.iterrows():
        migr[int(r1['node_label'])] = {}

        for i2, r2 in df.iterrows():
            if r2['node_label'] == r1['node_label']:
                pass
            else:
                d = vincenty((r1['lat'], r1['lon']), (r2['lat'], r2['lon'])).km
                migr[int(r1['node_label'])][int(r2['node_label'])] = compute_migr_prob(grav_params, r1['pop'], r2['pop'], d)

        p_sum[jj] = np.sum(migr[r1['node_label']].values())
        jj += 1

    return migr




grav_params = np.array([2.*7.50395776e-06,   9.65648371e-01,   9.65648371e-01, -1.10305489e+00])
grid_df = pd.read_csv("../../data/gridded_pop/cleaned/mapatizya_cells_jaline.csv")
migr_dict = compute_migr_dict(grid_df,grav_params)
with open("../../data/gridded_pop/cleaned/mapatizya_rates_jaline.json", 'w') as f:
    json.dump(migr_dict, f, indent=4)