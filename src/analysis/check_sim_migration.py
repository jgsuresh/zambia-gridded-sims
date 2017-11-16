# Compare migration within a certain region of a simulation to data and model migration of the same region:

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from geopy.distance import vincenty
from grid_ids_to_nodes import generate_latlong_lookup


def get_trips_within_region(trips_df,lat_bnds,long_bnds):
    def constrain_field(df,field,range):
        return np.logical_and(df[field] > range[0], df[field] < range[1])

    box1 = constrain_field(trips_df,'lat_start',lat_bnds)
    box2 = constrain_field(trips_df, 'lat_end', lat_bnds)
    box3 = constrain_field(trips_df,'long_start',long_bnds)
    box4 = constrain_field(trips_df, 'long_end', long_bnds)
    box = np.logical_and(np.logical_and(box1,box2),np.logical_and(box3,box4))

    return np.array(trips_df[box]['d'])

def get_trips_from_sim_reporter(sim_df):
    # Return array of distances of all trips in simulation

    foo = sim_df.groupby(['Time','Individual_ID'])
    foo = foo['Node_ID'].unique()
    foo = foo.reset_index()

    demo_fp_full = 'C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/experiments/gravity_test_v0/inputs/Demographics/MultiNode/demo.json'
    coords_lookup = generate_latlong_lookup(demo_fp_full)

    def calc_d(row):
        # Lookup the node IDs:
        ids = row[2]
        if len(ids) == 2:
            lat1,long1 = coords_lookup[ids[0]]
            lat2,long2 = coords_lookup[ids[1]]

            d = vincenty((lat1,long1),(lat2,long2)).km
            return d
        else:
            return -1

    d = foo.apply(lambda x: calc_d(x),axis=1)
    return d[d>=0]

def get_model_trip_predictions(grid_df, grav_params):
    # Return array of distances of all trips among grid cells, with weighting given by the gravity model
    lat = np.array(grid_df['lat'])
    long = np.array(grid_df['lon'])
    p = np.array(grid_df['pop'])

    n_pix = len(lat)
    lat1 = np.repeat(lat,n_pix)
    long1 = np.repeat(long, n_pix)
    p1 = np.repeat(p, n_pix)
    lat2 = np.tile(lat,n_pix)
    long2 = np.tile(long, n_pix)
    p2 = np.tile(p, n_pix)


    # Calculate centroid-ed distance:
    func = lambda x: vincenty(x[0:2],x[2:4]).km

    dat_columns = np.column_stack((lat1, long1, lat2, long2))
    d = np.array(map(func, dat_columns))

    w = grav_params[0] * p1**grav_params[1] * p2**grav_params[2] * d**grav_params[3]

    non_self = d != 0
    return [d[non_self],w[non_self]]



# def compute_all_possible_distance_combinations

def plot_trip_histograms(catch_name,sim_d,model_d,model_w,data_d,with_sim=False,bins=10,maxd=12):
    plt.figure()
    #
    # hist1,bin_edges1 = np.histogram(sim_d,bins=50)
    # hist2,bin_edges2 = np.histogram(model_d,weights=model_w,bins=50)

    if with_sim:
        plt.hist(sim_d,histtype='step',lw=2,label="Simulation",normed=True,range=(1,maxd),bins=bins,color='C0')
    plt.hist(model_d,weights=model_w,histtype='step', lw=2,ls='dashed', label="Gravity Model", normed=True,range=(1,maxd),bins=bins,color='C1')
    plt.hist(data_d[data_d > 1],histtype='step', lw=3, label="Data", normed=True,range=(1,maxd),bins=bins,color='C2')
    plt.legend()
    plt.xlabel("Distance (km)")
    plt.ylabel("Frequency (normed)")
    plt.title(catch_name)
    plt.savefig("../../results/migration/{}_trip_hist.pdf".format(catch_name))
    plt.show()

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    from matplotlib.colors import LinearSegmentedColormap
    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    # return base.from_list(cmap_name, color_list, N)
    return LinearSegmentedColormap.from_list(cmap_name, color_list, N)

def plot_grid_scatter(catch_name,grid_df):
    if catch_name == 'Chiyabi':
        s = 250
    elif catch_name == 'Mapatizya':
        s = 25
    elif catch_name =='all':
        s = 5
    else:
        s = 50
    plt.figure()
    cm = discrete_cmap(5,base_cmap=plt.cm.viridis)
    plt.scatter(grid_df['lon'],grid_df['lat'],marker='s',cmap=cm,c=np.log10(grid_df['pop']),s=s,edgecolors='black')
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.title(catch_name)
    cb = plt.colorbar()
    cb.set_label("Log10[Pop]")
    plt.savefig("../../results/migration/{}_grid_scatter.pdf".format(catch_name))
    plt.show()

def plot_total_trip_count_per_pixel(grid_df,grav_params):
    # Loop over each pixel.
    p = np.array(grid_df['pop'])
    lat = np.array(grid_df['lat'])
    long = np.array(grid_df['lon'])

    n_pix = len(p)
    ind = np.arange(n_pix)
    n_trips_day = np.zeros(n_pix)
    # Then loop over all other pixels.
    for i in np.arange(n_pix):
        other_pix = ind != i

        func = lambda x: vincenty((lat[i],long[i]),(x[0:2])).km
        dat_columns = np.column_stack((lat[other_pix],long[other_pix]))
        d = np.array(map(func,dat_columns))

        trips_per_link = grav_params[0] * p[i]**grav_params[1] * p[other_pix]**grav_params[2] * d**grav_params[3]

        n_trips_day[i] = np.sum(trips_per_link)

    plt.figure()
    plt.scatter(p,n_trips_day/p*4)
    plt.xlabel("Pixel Pop")
    plt.ylabel("(Daily trips to/from)/Pixel Pop")
    plt.show()


def run_all(catch_name,with_sim=False,maxd=12,bins=10):
    # Plot histogram of actual trip data, correctly normalized for coverage/round-linkage effects.  [hard]
    # Plot histogram of model trip data, given a set of pixels with a certain arrangement and populations.  [easy]
    # Plot histogram of corresponding DTK sim trip data, given the same set of nodes.  [medium]



    # Get histogram of trip distances from DTK simulation:
    if with_sim:
        sim_report = 'C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/experiments/gravity_test_v0/data/ReportEventRecorder.csv'
        sim_df = pd.read_csv(sim_report)
        sim_d = get_trips_from_sim_reporter(sim_df)
    else:
        sim_d = None

    # Get histogram of trip distances inferred from gravity model:
    grav_params = np.array([7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00])
    # grav_params = np.array([7.50395776e-06, 1., 1., -2.])
    grid_pop_csv_file = '../../data/gridded_pop/cleaned/{}_max_pop.csv'.format(catch_name.lower())
    grid_df = pd.read_csv(grid_pop_csv_file)
    model_d,model_w = get_model_trip_predictions(grid_df,grav_params)

    # Get histogram of trip distances directly from data in this region:
    trips_file = 'C:/Users/jsuresh/OneDrive - IDMOD/Code/zambia/migration/cluster_on_all_individuals/all_trips_rds1-10.csv'
    trips_df = pd.read_csv(trips_file)
    lat_bnds = [np.min(grid_df['lat']),np.max(grid_df['lat'])]
    long_bnds = [np.min(grid_df['lon']), np.max(grid_df['lon'])]
    data_d = get_trips_within_region(trips_df,lat_bnds,long_bnds)

    # Plot histograms side by side:
    plot_trip_histograms(catch_name.lower(),sim_d,model_d,model_w,data_d,with_sim=with_sim,maxd=maxd,bins=bins)

    plot_grid_scatter(catch_name,grid_df)



# if __name__=="main":
# run_all('all',with_sim=False,maxd=300,bins=50)
# run_all('Chiyabi',with_sim=True)
# run_all('Luumbo',with_sim=False,maxd=20,bins=15)
# run_all('Mapatizya',with_sim=False,maxd=20,bins=15)
# run_all('Munyumbwe',with_sim=False,maxd=20,bins=15)

# Get histogram of trip distances inferred from gravity model:
catch_name = "Mapatizya"
grav_params = np.array([7.50395776e-06, 9.65648371e-01, 9.65648371e-01, -1.10305489e+00])
# grav_params = np.array([7.50395776e-06, 1., 1., -2.])
grid_pop_csv_file = '../../data/gridded_pop/cleaned/{}_max_pop.csv'.format(catch_name.lower())
grid_df = pd.read_csv(grid_pop_csv_file)

plot_total_trip_count_per_pixel(grid_df,grav_params)