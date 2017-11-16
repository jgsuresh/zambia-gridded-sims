import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
import seaborn
from struct import pack
import re
import dtk.tools.demographics.compiledemog as compiledemog
# from Mig

# Script which takes demographics file as input, and returns migration json file, using a gravity model fit.
from geopy.distance import vincenty
from sklearn.cluster import DBSCAN


def load_demo(demo_file):
    with open(demo_file, 'r') as f:
            demo_dict = json.load(f)

    N = demo_dict['Metadata']['NodeCount']
    lat = np.ones(N)
    long = np.ones(N)
    grid_id = np.ones(N)
    node_id = np.ones(N)
    pop = np.ones(N)

    for i in xrange(N):
        node = demo_dict['Nodes'][i]
        lat[i] = node['NodeAttributes']['Latitude']
        long[i] = node['NodeAttributes']['Longitude']
        grid_id[i] = node['NodeAttributes']['FacilityName']
        node_id[i] = node['NodeID']
        pop[i] = node['NodeAttributes']['InitialPopulation']

    df = pd.DataFrame({
        'lat': lat,
        'long': long,
        'grid_id': grid_id,
        'node_id': node_id,
        'pop': pop
    })

    return df

# def cluster_pops_v2(df,plot=True):
    # Distribute population of each pixel randomly within the pixel, then do a DBSCAN clustering
    # def distribute_pop:
    # def get_width()
    # dlat =
    # dlong =

    # return labels

def cluster_pops(df, plot=True):
    print "DBSCAN clustering..."
    clusterer = DBSCAN(min_samples=5, eps=0.0095)
    clusterer.fit(np.column_stack((df['lat'],df['long'])),sample_weight=df['pop'])
    labels = np.copy(clusterer.labels_)
    print "DBSCAN clustering finished."

    if plot:
        plt.figure()
        for lb in np.unique(labels):
            in_lb = labels==lb
            color = None
            if lb == -1:
                color = 'black'
            plt.scatter(df['long'][in_lb],df['lat'][in_lb],c=color,s=df['pop'][in_lb])
        plt.show()
    return labels

def compute_migr_prob(grav_params,ph,pd,d):
    num_trips = grav_params[0] * ph**grav_params[1] * pd**grav_params[2] * d**grav_params[3]
    prob_trip = np.min([1.,num_trips/ph])
    return prob_trip

def compute_migr_dict(df, grav_params, d_thresh=-1, return_prob_sums=False):
    migr = {}

    p_sum = np.zeros(len(df))
    jj = 0
    for i1,r1 in df.iterrows():
        migr[r1['node_id']] = {}

        for i2,r2 in df.iterrows():
            if r2['node_id'] == r1['node_id']:
                pass
            else:
                d = vincenty((r1['lat'],r1['long']),(r2['lat'],r2['long'])).km
                migr[r1['node_id']][r2['node_id']] = compute_migr_prob(grav_params,r1['pop'],r2['pop'],d)

        p_sum[jj] = np.sum(migr[r1['node_id']].values())
        jj += 1

    if return_prob_sums:
        return [migr,p_sum]
    else:
        return migr

'''
save link rates to a human readable file;
the txt file is consumable by link_rates_txt_2_bin(self) like function to generate DTK migration binary
'''


def save_link_rates_to_txt(rates_txt_file_path, link_rates):
    with open(rates_txt_file_path, 'w') as fout:
        for src, v in link_rates.items():
            for dest, mig in v.items():
                fout.write('%d %d %0.1g\n' % (int(src), int(dest), mig))

'''
convert a txt links rates file (e.g. as generated by save_link_rates_to_txt(self...)) to DTK binary migration file 
'''

def link_rates_txt_2_bin(rates_txt_file_path, rates_bin_file_path, route="local"):

    fopen = open(rates_txt_file_path)
    fout = open(rates_bin_file_path, 'wb')

    net = {}
    net_rate = {}

    MAX_DESTINATIONS_BY_ROUTE = {'local': 100,
                                 'regional': 30,
                                 'sea': 5,
                                 'air': 60}

    for line in fopen:
        s = line.strip().split()
        ID1 = int(float(s[0]))
        ID2 = int(float(s[1]))
        rate = float(s[2])
        # print(ID1,ID2,rate)
        if ID1 not in net:
            net[ID1] = []
            net_rate[ID1] = []
        net[ID1].append(ID2)
        net_rate[ID1].append(rate)

    for ID in sorted(net.keys()):

        ID_write = []
        ID_rate_write = []

        if len(net[ID]) > MAX_DESTINATIONS_BY_ROUTE[route]:
            print('There are %d destinations from ID=%d.  Trimming to %d (%s migration max) with largest rates.' % (
            len(net[ID]), ID, MAX_DESTINATIONS_BY_ROUTE[route], route))
            dest_rates = zip(net[ID], net_rate[ID])
            dest_rates.sort(key=lambda tup: tup[1], reverse=True)
            trimmed_rates = dest_rates[:MAX_DESTINATIONS_BY_ROUTE[route]]
            # print(len(trimmed_rates))
            (net[ID], net_rate[ID]) = zip(*trimmed_rates)
            # print(net[ID], net_rate[ID])

        for i in xrange(MAX_DESTINATIONS_BY_ROUTE[route]):
            ID_write.append(0)
            ID_rate_write.append(0)
        for i in xrange(len(net[ID])):
            ID_write[i] = net[ID][i]
            ID_rate_write[i] = net_rate[ID][i]
        s_write = pack('L' * len(ID_write), *ID_write)
        s_rate_write = pack('d' * len(ID_rate_write), *ID_rate_write)
        fout.write(s_write)
        fout.write(s_rate_write)

    fopen.close()
    fout.close()

def save_migration_header(demographics_file_path, outfilename=None):

    # generate migration header for DTK consumption
    # todo: the script below needs to be refactored/rewritten
    # in its current form it requires compiled demographisc file (that's not the only problem with its design)
    # to compile the demographics file need to know about compiledemog file here, which is unnecessary
    # compiledemog.py too could be refactored towards object-orientedness
    # the demographics_file_path supplied here may be different from self.demographics_file_path)
    compiledemog.main(demographics_file_path)
    import createmigrationheader
    createmigrationheader.main('dtk-tools', re.sub('\.json$', '.compiled.json', demographics_file_path), 'local',
                               outfilename=outfilename)


def gen_gravity_links_json(demo_file,grav_params,outf=None):
    df = load_demo(demo_file)
    migr_dict = compute_migr_dict(df, grav_params, return_prob_sums=False)

    # Save to file:
    if outf==None:
        outf = 'grav_migr_rates.json'
    with open(outf, 'w') as f:
        json.dump(migr_dict, f, indent=4)

    return migr_dict



if __name__=="__main__":
    grid_pop_csv_file = 'C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/gridding/test-10-17/grid_lookup_max_pop_refmt.csv'
    demo_file = 'C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/experiments/gravity_test_v0/inputs/Demographics/MultiNode/test_demographics.json'

    df = load_demo(demo_file)
    # labels = cluster_pops(df)

    # Treat each pixel as an independent cluster, since gravity model is roughly linear in population
    # grav_params = np.array([4.60248607e-05,   9.04960570e-01,   9.04960523e-01,-1.28511700e+00])
    grav_params = np.array([7.50395776e-06,   9.65648371e-01,   9.65648371e-01, -1.10305489e+00])

    migr_dict,p_sum = compute_migr_dict(df, grav_params, return_prob_sums=True)

    # Save to file:
    with open('migr_test.json', 'w') as f:
        json.dump(migr_dict, f, indent=4)

    save_link_rates_to_txt('test.txt',migr_dict)

    link_rates_txt_2_bin('test.txt','test_migration.bin')

    save_migration_header(demo_file)