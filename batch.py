## IMPORTS

# Numerical packages and methods
import numpy as np
from sklearn import linear_model
from scipy.optimize import curve_fit
import scipy.stats as st
import scipy.interpolate as interp
import pandas as pd
import statsmodels.nonparametric.smoothers_lowess as lowess

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
import sys
import progressbar



# Import what ?
prefix = "./data/" + sys.argv[1] + "/"
directory = [f for f in os.listdir(prefix) if f.endswith(".csv")]
names = [i.split(".")[0].capitalize() for i in directory]


# Params
sensitivity = 0
periodicity = 24
penalty = 1e-3
delay = 8
vaccine = 1965


# Results directory
if not os.path.isdir(prefix + "results") :
	os.mkdir(prefix + "results")














# Now, for each file in the directory
for idx, file in enumerate(directory) :



	# Import data
	data = pd.read_csv(prefix + file)

	# Time and vaccinations
	t = data["time"].values
	v = np.where(t > 1965)[0][0] if t[-1] > vaccine else len(t)
	t = t[:v]

	# Births and cases
	B = data["births"].values[:v]
	C = data["reported_cases"].values[:v]






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
	plt.figure()

	plt.subplot(211)
	plt.plot(t, C, linewidth=2)
	plt.title("Reported Cases, %s" % names[idx])
	plt.xlabel("Time")
	plt.ylabel("Cases")

	plt.subplot(212)
	plt.plot(t, B, linewidth=2)
	plt.title("Live Births, %s" % names[idx])
	plt.xlabel("Time")
	plt.ylabel("Births")

	plt.tight_layout()
	plt.savefig(prefix + "results/%s_0_timeseries.pdf" % names[idx])
	print "%s Time-Series done." % names[idx]
















	# Susceptible Reconstruction

	# Compute cumulative births and incidence
	Y = np.cumsum(B)
	X = np.cumsum(C) 



	"""
	# Rho
	Yhat = lowess.lowess(Y, X, 0.7)[:, 1]

	# Now sample it with indices
	yhat = interp.interp1d(X, Yhat)

	# Correct the first value
	x = np.linspace(0, X[-1], len(X))
	yhat = yhat(x)
	yhat[0] = 2 * yhat[1] - yhat[2]

	# Rho : put some splines through Yhat and take the derivative
	rho = interp.UnivariateSpline(x, yhat).derivative()(x)
	
	"""
	reg = linear_model.BayesianRidge(fit_intercept=False, compute_score=True)

	# Compute the R^2 for a range of polynomials from degree-1 to degree-7
	# The fit score has a penalty proportional to the square of the degree of the polynomial
	"""
	Ns = range(2, 8)
	#scores = []
	for n in Ns :
	    reg.fit(np.vander(X, n), Y)
	    scores.append(reg.score(np.vander(X, n), Y) - penalty * n**2)
	"""    
	# Use the polynomial that maximised R^2 to compute Yhat
	Yhat = reg.fit(np.vander(X, 4), Y).predict(np.vander(X, 4))
	
	
	# Compute rho as the derivative of the splines that are fit between X and the estimated Y
	rho = interp.UnivariateSpline(X, Yhat).derivative()(X)
	
	# Compute Z as the residuals of regression
	Z = Y - Yhat





	# Plots
	plt.figure()

	plt.subplot(221)
	plt.plot(t, X, linewidth=2)
	plt.plot(t, Y, linewidth=2)
	plt.plot(t, Yhat, linewidth=2)
	plt.title("Reported and Inferred Cases, %s" % names[idx])
	plt.legend(["Reported Cases", "Cumulative Births", "Inferred Cases"], loc=2)

	plt.subplot(222)
	plt.axhline(1./np.mean(rho), color="r", linewidth=2)
	plt.plot(t, 1./rho, linewidth=2)
	plt.ylim([0, np.max(1.1/rho)])
	plt.title(r"Inferred Reporting Rate $1/\rho_t = %.2f$" % (1./np.mean(rho)))
	#plt.legend([r"$E[1/\rho_t]=%.2f$" % (1./np.mean(rho))])

	plt.subplot2grid((2, 2), (1, 0), colspan=2)
	plt.plot(t, Z, linewidth=2)
	plt.title("Susceptible Dynamics $Z_t$")
	plt.xlabel("Time (years)")

	"""
	plt.plot(np.array(Ns)-1, scores, linewidth=2)
	plt.axvline(Ns[np.argmax(scores)]-1, color="r", linewidth=2)
	plt.title("Polynomial Model Fit, n = %d" % (Ns[np.argmax(scores)]-1))
	plt.xlabel("Polynomial Degree")
	plt.ylabel("Penalised Goodness of Fit")
	"""
	plt.tight_layout()
	plt.savefig(prefix + "results/%s_1_susceptible_reconstruction.pdf" % names[idx])
	print "%s Susceptible Reconstruction done." % names[idx]
	















	# Fit Sbar

	# All possible values of Sbar
	Svals = np.linspace(np.abs(np.min(Z))+1, np.abs(np.min(Z))*20, 200)

	# Likelihood of fit
	l = np.zeros(len(Svals))



	# Define our parameters
	params = lmfit.Parameters()
	params.add("alpha", min=0.5, max=.9999, value=0.95) # Alpha
	for i in range(periodicity) : # Seasonalities
	    params.add("r%d" % i, value=0.)
	rstr = ["r%d" % (i % periodicity) for i in list(itertools.chain.from_iterable(epi))][:-1]
	#if 
	    
	# Objective function
	def profile_residuals(params, rho, C, Z, z, Sestimate) :
	    alphafit = params["alpha"].value
	    r = [params[i].value for i in rstr]
	    if np.isnan(Sestimate) :
	        Sestimate = params["Sest"].value
	    return alphafit * np.log(rho[z[:-1]]*C[z[:-1]]) + r + np.log(Sestimate + Z[z[:-1]]) - np.log(rho[z[1:]]*C[z[1:]])
	    

	    
	    
	    
	    
	    
	# Compute best fit for each possible Sbar
	for i, Sestimate in enumerate(Svals) :
	    l[i] = lmfit.minimize(profile_residuals, params, args=(rho, C, Z, z, Sestimate), method="leastsq").chisqr
	    
	    
	# Fit window
	fitwindow = 15
	fitwindowL = np.min([fitwindow, np.argmin(l)])
	fitwindowR = np.min([fitwindow, len(Svals) - np.argmin(l)])
	    
	# Run again using scan estimate
	params.add("Sest", value = Svals[np.argmin(l)])
	L = lmfit.minimize(profile_residuals, params, args=(rho, C, Z, z, np.nan), method="leastsq")


	# Extract parameters and errors
	Sbar = L.params["Sest"].value
	r = np.exp([L.params["r" + str(i)].value for i in range(periodicity)])
	alphaSbar = L.params["alpha"].value
	errup = np.exp(np.log(r) + [2*L.params["r" + str(i)].stderr for i in range(periodicity)])
	errdn = np.exp(np.log(r) - [2*L.params["r" + str(i)].stderr for i in range(periodicity)])
	    
	    



	# Plot
	plt.figure()

	plt.subplot(121)
	plt.axvline(x=Sbar, color="red", linewidth=2)
	plt.loglog(Svals, l, linewidth=2)
	plt.title("Goodness of Fit, %s" % names[idx])
	plt.xlabel(r"$\bar{S}$")
	plt.ylabel(r"$\chi^2$")
	plt.legend([r"Profile Likelihood $\bar{S}$ = %d" % Sbar])


	plt.subplot(122)
	plt.plot(r, linewidth=3)
	plt.fill_between(range(periodicity), errup, errdn, color=colours[0], alpha=0.3)
	plt.xlim([0, periodicity])
	plt.title("Periodicity")
	plt.xlabel("Period")
	
	plt.tight_layout()
	plt.savefig(prefix + "results/%s_2_meanS_periodicity.pdf" % names[idx])
	print "%s Mean S and Periodicity done." % names[idx]














	# Simulations

	allss = []
	allsd = []
	allrs = []
	allrd = []
	I = []
	pbar = progressbar.ProgressBar(widgets = \
			[progressbar.Percentage(), progressbar.Bar(), progressbar.ETA()], \
			maxval = int(sys.argv[2]))
	pbar.start()

	for q in range(int(sys.argv[2])) :
	    
	    predI = np.zeros_like(C)
	    predS = np.zeros_like(C)
	    starts = [e[0] for e in epi]
	    
	    # Epi size and duration
	    simsizes = []
	    simduration = []
	    realsizes = []
	    realduration = []
	    maxdur = np.append(np.diff(starts), np.inf)
	    
	    # Seed initial epidemic points
	    for index, e in enumerate(starts) :
	        predI[e] = np.ceil(rho[e] * C[e])
	        predS[e] = np.ceil(Sbar + Z[e])
	    
	        for i in range(e+1, len(C)) :
	            predI[i] = np.random.negative_binomial(max(np.round(predI[i-1]), 1), \
	                                                   max(np.round(predI[i-1]), 1) / ( max(np.round(predI[i-1]), 1) + \
	                                                               r[i % periodicity] * ( predI[i-1] ** alphaSbar ) * predS[i-1]))
	            #predI[i] = np.random.poisson(r[i % periodicity] * ( predI[i-1] ** alphaSbar ) * predS[i-1])
	            predS[i] = B[max(i - delay, 0)] + predS[i-1] - predI[i]
	            
	        simsizes.append(np.sum(predI[e:]))
	        simduration.append(min(np.argmin(predI[e:]), maxdur[index]))
	        realduration.append(np.argmin(C[e:] * rho[e:]))
	        realsizes.append(np.sum(rho[e:e+realduration[-1]] * C[e:e+realduration[-1]]))
	        
	        pbar.update(q)

	    allss.append(simsizes)
	    allrs.append(realsizes)
	    allsd.append(simduration)
	    allrd.append(realduration)
	    I.append(predI)
	    

	    
	s0plotx2 = []
	s0ploty2 = []
	for j in allss :
	    for i, e in enumerate(starts) :
	        if j[i] > 250 :
	            s0plotx2.append(Sbar + Z[e])
	            s0ploty2.append(j[i])
	    
	    
	newshape = np.array(allss).shape[0] * np.array(allss).shape[1]
	allss = np.array(allss).ravel().reshape(newshape, 1)
	allrs = np.array(allrs).ravel().reshape(newshape, 1)
	allsd = np.array(allsd).ravel().reshape(newshape, 1)
	allrd = np.array(allrd).ravel().reshape(newshape, 1)
	    
	#idx = allrs > 500
	#allss = allss[idx].reshape(np.sum(idx), 1)
	#allrs = allrs[idx].reshape(np.sum(idx), 1)

	    
	sfit = linear_model.BayesianRidge(fit_intercept=False)
	dfit = linear_model.BayesianRidge(fit_intercept=False)

	dfit.fit(allrd, allsd)
	sfit.fit(allrs, allss)

	I = np.array(I)    
	    
	    
	# Plot   
	
	plt.figure()
	for i in I :
	    plt.plot(i, c = colours[2], linewidth=1, alpha=0.01)
	#plt.plot(np.mean(I, axis=0), c=colours[2], linewidth=2)
	plt.plot(C*rho, c = colours[0], linewidth=4, alpha=0.6)

	plt.tight_layout()
	plt.savefig(prefix + "results/%s_3_predictions.pdf" % names[idx])
	print "%s Predictions done." % names[idx]





















	# Size and Duration of Epidemics : Real vs Predicted
	plt.figure()

	plt.subplot(211)
	plt.title("%s, Size of Epidemics : Gradient = %f, R^2 = %f" % (names[idx], sfit.coef_[0], sfit.score(np.array(realsizes).reshape(len(realsizes), 1), realsizes)))
	plt.xlabel("Real Size")
	plt.ylabel("Simulated Size")
	for i in range(q) :
	    plt.scatter(allrs[i], allss[i], alpha=0.3, c=colours[2], s=35)
	plt.plot(realsizes, sfit.predict(np.array(realsizes).reshape(len(realsizes), 1)), linewidth=2)
	    
	plt.subplot(212)
	plt.title("Duration of Epidemics, Gradient = %f, R^2 = %f" % (dfit.coef_[0], dfit.score(np.array(realduration).reshape(len(realduration), 1), realduration)))
	plt.xlabel("Real Duration")
	plt.ylabel("Simulated Duration")
	for i in range(q) :
	    plt.scatter(allrd[i], allsd[i], alpha=0.3, c=colours[2], s=35)
	plt.plot(realduration, dfit.predict(np.array(realduration).reshape(len(realduration), 1)), linewidth=2)

	plt.tight_layout()
	plt.savefig(prefix + "results/%s_4_sizes_durations.pdf" % names[idx])
	print "%s Sizes and Durations done." % names[idx]





















	# Susceptibles vs Sizes 

	s0plotx = []
	s0ploty = []
	for i, e in enumerate(starts) :
	    if realsizes[i] > 250 :
	        s0plotx.append(Sbar + Z[e])
	        s0ploty.append(realsizes[i])
	        
	s0 = linear_model.BayesianRidge(fit_intercept=True)
	s0.fit(np.array(s0plotx).reshape(len(s0plotx), 1), np.array(s0ploty).reshape(len(s0plotx), 1))
	            
	s1 = linear_model.BayesianRidge()
	s1.fit(np.array(s0plotx2).reshape(len(s0plotx2), 1), np.array(s0ploty2).reshape(len(s0plotx2), 1))



	# Plot
	plt.figure()

	plt.subplot(211)
	plt.scatter(s0plotx, s0ploty, c = colours[0])
	plt.plot(s0plotx, s0.predict(np.array(s0plotx).reshape(len(s0plotx), 1)), linewidth=2)
	plt.title("%s S0 vs Real Epidemic Size, Gradient = %f, Intercept = %f" % (names[idx], s0.coef_[0], s0.predict(0)))


	plt.subplot(212)
	plt.scatter(s0plotx2, s0ploty2, c = colours[0], alpha=0.3)
	plt.plot(s0plotx2, s1.predict(np.array(s0plotx2).reshape(len(s0plotx2), 1)), linewidth=2)
	plt.title("S0 vs Simulated Epidemic Size, Gradient = %f, Intercept = %f" % (s1.coef_[0], s1.predict(0)))

	plt.tight_layout()
	plt.savefig(prefix + "results/%s_5_s0_vs_size.pdf" % names[idx])
	print "%s S0 vs Sizes done." % names[idx]














