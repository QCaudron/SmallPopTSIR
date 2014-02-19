import numpy as np
import scipy.interpolate as interp
import pandas as pd
import os


# Create directory if needed
if not os.path.isdir("../data/uk_all") :
	os.mkdir("../data/uk_all")


# Import raw data
pop = pd.read_csv("../raw_data/uk/population_uk.csv")
births = pd.read_csv("../raw_data/uk/births_uk.csv")
cases = pd.read_csv("../raw_data/uk/cases_uk.csv")


# Determine columns
columns = pop.columns[1:]


# Time array - from previously-calculated London data, to ensure we stay constant
time = pd.read_csv("../data/uk_big/london.csv")["time"].values


# Iterate over columns
for i in columns :

	# Population
	p = interp.interp1d(pop.values[:, 0], pop[i].values, kind = "cubic")(time)

	# Births
	b = interp.interp1d(births.values[:, 0], births[i].values, kind = "cubic")(time)
	b /= np.sum(b)
	b *= np.sum(births[i].values)

	# Reported Cases
	c = interp.interp1d(cases.values[:, 0], cases[i].values, kind = "linear")(time)
	c /= np.sum(c)
	c *= np.sum(cases[i].values[:np.where(cases.values[:, 0] > time[-1])[0][0]])

	# Generate output dataframe
	out = {}
	out["births"] = np.round(b).astype(int)
	out["population"] = np.round(p).astype(int)
	out["reported_cases"] = np.round(c).astype(int)
	out["time"] = time

	pd.DataFrame(out).to_csv("../data/uk_all/" + i.lower() + ".csv", \
		 header = True, cols = ["time", "births", "population", "reported_cases"], index = False)

	# Print the name of the town
	print i.capitalize()