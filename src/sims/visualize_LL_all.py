import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid")

base = 'C:/Users/jsuresh/OneDrive - IDMOD/Projects/zambia-gridded-sims/'
LL_all_path = base + "src/sims/Serial_Motaze/_plots/LL_all.csv"

LL_all = pd.read_csv(LL_all_path)
plt.scatter(LL_all['arabiensis_scale'],LL_all['arabiensis_funestus_ratio'],c=LL_all['total'])
plt.xlabel('arabiensis_scale')
plt.ylabel('arabiensis_funestus_ratio')
plt.colorbar()
plt.show()