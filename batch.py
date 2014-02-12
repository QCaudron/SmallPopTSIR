## IMPORTS

# Numerical packages and methods
import numpy as np
from sklearn import linear_model
from scipy.optimize import curve_fit
import scipy.stats as st
import scipy.interpolate as interp
import pandas as pd

# Monte Carlo and Nonlinear Fitting
import lmfit

# Plotting
import matplotlib.pyplot as plt
import matplotlib
import seaborn
colours = seaborn.color_palette("deep", 8)

# Other
import itertools
import os




# Import what ?
prefix = "./data/" + argv[1] + "/"
directory = os.listdir(prefix)
names = [i.split(".")[0] for i in directory]


# Params
sensitivity = 0
periodicity = 24


# Results directory
if not os.path.isdir(prefix + "results") :
	os.mkdir(prefix + "results")
	



# Now, for each file in the directory
for idx, file in enumerate(directory) :



	# Import data
	data = pd.read_csv(prefix + file)

	# Time, births and cases
	t = data["time"].values
	B = data["births"].values
	C = data["cases"].values






	# Find epidemics
	epi = []

	# If there are more than 20% zeros
	if (np.sum(C <= sensitivity).astype(float) / len(C)) > 0.2 :
	    z = np.where(C > sensitivity)[0] # Find epidemics over sensitivity threshold
	    dz = np.where(np.append(np.insert(np.diff(z), 0, 0), -1) != 1)[0]
	    for i in range(len(dz)-1) :
	        epi.append(z[dz[i]:dz[i+1]])
	      
	else : # Otherwise, slice at local minima using smoothed zero-crossings in the derivative
	    z = range(len(C))
	    z2 = np.diff(np.convolve(C, np.hanning(19), "same"))
	    dz = np.append(np.insert((np.where((z2[:-1] < 0) * (z2[1:] > 0) == True)[0]), 0, 0), len(C))
	    for i in range(len(dz)-1) :
	        epi.append(range(dz[i], dz[i+1]))

	epi = np.array(epi)






	# Plot
	f = plt.figure()

	plt.subplot(211)
	plt.plot(t, C, linewidth=3)
	plt.title("Reported Cases, %s" % names[idx])
	plt.xlabel("Time")
	plt.yabel("Cases")

	plt.subplot(212)
	plt.plot(t, B, linewidth=3)
	plt.title("Live Births, %s" % names[idx])
	plt.xlabel("Time")
	plt.yabel("Births")

	plt.tight_layout()
	plt.savefig(prefix + "results/")







	#############################################
	# Susceptible Reconstruction

	# Compute cumulative births and incidence
	Y = np.cumsum(B)
	X = np.cumsum(C) 

	# Rho
	reg = linear_model.BayesianRidge(fit_intercept=False, compute_score=True)

	# Compute the R^2 for a range of polynomials from degree-1 to degree-10
	# The fit score has a penalty proportional to the square of the degree of the polynomial
	Ns = range(2, 12)
	scores = []
	for n in Ns :
	    reg.fit(np.vander(X, n), Y)
	    scores.append(reg.score(np.vander(X, n), Y) - penalty * n**2)
	    
	# Use the polynomial that maximised R^2 to compute Yhat
	Yhat = reg.fit(np.vander(X, Ns[np.argmax(scores)]), Y).predict(np.vander(X, Ns[np.argmax(scores)]))

	# Compute rho as the derivative of the splines that are fit between X and the estimated Y
	rho = interp.UnivariateSpline(X, Yhat).derivative()(X)

	# Compute Z as the residuals of regression
	Z = Y - Yhat





























