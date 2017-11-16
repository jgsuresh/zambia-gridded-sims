import pandas as pd

def convert_cait_CSV(cait_CSV, outname):
    new_CSV = cait_CSV.copy()
    new_CSV = new_CSV.rename(columns={'maxpop': 'pop','grid.cell':'node_label','mid.x': 'lon', 'mid.y': 'lat'})
    new_CSV = new_CSV.drop('catch',axis=1)
    new_CSV.to_csv(outname)

if __name__=="__main__":
    df = pd.read_csv('grid_lookup_max_pop.csv')
    convert_cait_CSV(df,'grid_lookup_max_pop_refmt.csv')