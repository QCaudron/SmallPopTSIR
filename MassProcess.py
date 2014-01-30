# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

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
matplotlib.rcParams['savefig.dpi'] = 1.5 * matplotlib.rcParams['savefig.dpi']
figsize(14, 8)
import seaborn
colours = seaborn.color_palette("deep", 8)

# Other
import itertools
import os

# <codecell>

# Import what ?
UK = 0
ICELAND = 1
FAROE = 0
BORNHOLM = 0


# Show graphs, or dump to file ?
visual = 1

# <codecell>

# Model parameters
periodicity = 24
delay = 8
sensitivity = 0
penalty = 5e-4




# Which data to import ?

if ICELAND :
    prefix = "./data/iceland/"
    directory = os.listdir(prefix)
elif LONDON :
    prefix = "./data/uk/"
    directory = os.listdir("./data/uk/")
elif FAROE :
    prefix = "./data/faroe/"
    directory = os.listdir("./data/faroe/")
elif BORNHOLM :
    prefix = "./data/bornholm/"
    directory = os.listdir("./data/bornholm/")



# Import births and cases

names = [i.split(".")[0] for i in directory]
data = []
B = []
C = []

for i in directory :
    data.append(pd.read_csv(prefix + i))
    B.append(data
        
        
        

        
    # Careful interpolation of spiky data
    
    # First, trim the data vector to keep an integer number of years
    data = np.delete(data, np.s_[np.where(data[:, 1] == np.floor(data[-1, 1]))[0] :], 0)
    
    # Define desired times : <periodicity> per year, to maintain 12-multiple exactly
    t = np.linspace(data[0, 1], data[-1, 1], np.ceil(data[-1, 1] - data[0, 1]) * periodicity - 1)
    
    
    
    
    # Births (B) and Cases (C) use simple linear interpolation, except London, which is already biweekly
    B = np.interp(t, data[:, 1], data[:, 2])
    if not LONDON :
        B = B / periodicity # births were reported annually
    C = np.interp(t, data[:, 1], data[:, 0]) 
    if not LONDON :
        C = np.round(C / periodicity * 12) # cases were reported monthly; round for integer cases
        
      


        
# Where are the epidemics ?
epi = []

# If there are many zeros ( here, we say at least 30% ), we can cut epidemics naturally at non-zeros
if (np.sum(C <= sensitivity).astype(float) / len(C)) > 0.3 :
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
        
        
        
# Plots
subplot(211)
plt.plot(t, C, linewidth=3)
#for e in epi :
#    axvline(t[e[0]], color="red", linewidth=2)
#    axvline(t[e[-1]], color="red", linewidth=2)
#    axhline(1.04*np.max(C), xmin=(t[e[0]]-t[0])/(t[-1]-t[0]), xmax=(t[e[-1]]-t[0])/(t[-1]-t[0]), linewidth=3, color="red")
title("Observed Cases, $C_t$. Red lines delimit epidemic.")
xlim([t[0], t[-1]])
ylim([-5, 1.05*np.max(C)])
xlabel("Time (index)")
ylabel("Cases")

subplot(212)
title("Birth Rates, $B_i$")
xlabel("Time (Years)")
ylabel("Births")
plt.plot(t, B, linewidth=3)

tight_layout()

# <codecell>

np.sum(C)

# <codecell>

# Susceptible Reconstruction

# Compute cumulative births and incidence
Y = np.cumsum(B)
X = np.cumsum(C) 


# Compute rho ( rate of reporting ) using Bayesian ridge regression with a polynomial
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



# Plots
subplot(221)
plt.plot(t, X, linewidth=3)
plt.plot(t, Y, linewidth=3)
plt.plot(t, Yhat, linewidth=3)
title("Reported and Inferred Cases")
legend(["Reported Cases", "Cumulative Births", "Inferred Cases"], loc=2)

subplot(222)
axhline(1./np.mean(rho), color=colours[2], linewidth=3)
plt.plot(t, 1./rho, linewidth=3)
ylim([0, 1])
title(r"Inferred Reporting Rate $1/\rho_t$")
legend([r"$E[1/\rho_t]=$" + str(1./np.mean(rho))])

subplot(223)
plt.plot(t, Z, linewidth=3)
title("Susceptible Dynamics $Z_t$")
xlabel("Time (years)")

subplot(224)
plt.plot(np.array(Ns)-1, scores, linewidth=3)
axvline(Ns[np.argmax(scores)]-1, color=colours[2], linewidth=3)
title("Polynomial Model Fit, n = " + str(Ns[np.argmax(scores)]-1))
xlabel("Polynomial Degree")
ylabel("Penalised Goodness of Fit")

plt.tight_layout()

# <codecell>

# EQUATION 15

# Fit a linear model to infer periodicity, alpha, and Sbar - using Z only

# Allocate design matrix
A = np.zeros((len(z)-1, periodicity+2))

# Periodicity indicators for the design matrix
for i in range(len(z)-1) :
    A[i, i % periodicity] = 1

# Set I(t-1), Z(t-1)
A[:, periodicity] = np.log(np.ceil(rho[z[:-1]] * C[z[:-1]]))
A[:, periodicity+1] = Z[z[:-1]]



# Initialise results vector
y = np.log(np.ceil(rho[z[1:]] * C[z[1:]]))



# Infer parameters using Bayesian ridge regression
reg2 = linear_model.BayesianRidge(fit_intercept=False)
reg2.fit(A, y)


# Extract useful parameters
rstar = np.exp(reg2.coef_[:periodicity]) # Sbar * r_t
alphaZ = reg2.coef_[periodicity] # alpha
zeta = reg2.coef_[periodicity+1] # Sbar

#plt.plot(rstar, linewidth=2)
#title("Periodicity")
print "Alpha = " + str(alphaZ)
print "Sbar = " + str(1./zeta)
if not (LONDON or ICELAND or FAROE or BORNHOLM) :
    print "Real Sbar = " + str(np.mean(S))
    print "Error in Sbar prediction = " + str(100*np.abs(1./zeta - np.mean(S))/np.mean(S)) + str(" %")

# <codecell>

# EQUATION 12

# All possible values of Sbar
Svals = np.linspace(np.abs(np.min(Z))+1, np.abs(np.min(Z))*13, 100)


# Likelihood of fit
l = np.zeros(len(Svals))


# Define our parameters
params = lmfit.Parameters()
params.add("alpha", min=0.5, max=1., value=0.95) # Alpha
for i in range(periodicity) : # Seasonalities
    params.add("r" + str(i), value=0.)
rstr = ["r" + str(i % periodicity) for i in list(itertools.chain.from_iterable(epi))][:-1]
#if 
    
# Objective function
def profile_residuals(params, rho, C, Z, z, Sestimate) :
    alphafit = params["alpha"].value
    r = [params[i].value for i in rstr]
    if isnan(Sestimate) :
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
subplot(121)
plt.axvline(x=Sbar, c=colours[1], linewidth=3)
plt.axvline(x=1/zeta, c=colours[2], linewidth=3)
plt.loglog(Svals, l, linewidth=3)
title("Goodness of Fit")
xlabel(r"$\bar{S}$")
ylabel(r"$\chi^2$")
legend([r"Profile Likelihood $\bar{S}$ = " + str(int(Sbar)), r"Taylor Expansion $\bar{S}$ = " + str(int(1/zeta))])


subplot(122)
#plt.plot(rstar, linewidth=2)
plt.plot(r, linewidth=3)
plt.fill_between(range(periodicity), errup, errdn, color=colours[2], alpha=0.3)
plt.plot(rstar * zeta, linewidth=3)
xlim([0, periodicity])
title("Periodicity")
xlabel("Period")
legend(["Using Full Equation", "Using Taylor Expansion"])
tight_layout()

# <codecell>

allss = []
allsd = []
allrs = []
allrd = []
I = np.zeros_like(C)


for q in range(1000) :
    
    # Initialise
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
        
        
    allss.append(simsizes)
    allrs.append(realsizes)
    allsd.append(simduration)
    allrd.append(realduration)
    I += predI
    

    
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

    
    
    
# Plot   
plt.subplot(211)
plt.plot(C*rho, linewidth=3)
plt.plot(I/q, c=colours[2], linewidth=3)
title("Per-Epidemic Predictions")
legend(["Observed", "Predicted"], loc=2)

plt.subplot(212)
plt.plot(predI - C*rho, linewidth=3)
title("Errors")

plt.tight_layout()

# <codecell>

plt.figure()
subplot(211)
plt.title("Size of Epidemics, Gradient = %f, R^2 = %f" % (sfit.coef_[0], sfit.score(np.array(realsizes).reshape(len(realsizes), 1), realsizes)))
plt.xlabel("Real Size")
plt.ylabel("Simulated Size")
for i in range(q) :
    plt.scatter(allrs[i], allss[i], alpha=0.3, c=seaborn.color_palette("deep", 3)[2], s=35)
plt.plot(realsizes, sfit.predict(np.array(realsizes).reshape(len(realsizes), 1)), linewidth=3)
    
subplot(212)
plt.title("Duration of Epidemics, Gradient = %f, R^2 = %f" % (dfit.coef_[0], dfit.score(np.array(realduration).reshape(len(realduration), 1), realduration)))
plt.xlabel("Real Duration")
plt.ylabel("Simulated Duration")
for i in range(q) :
    plt.scatter(allrd[i], allsd[i], alpha=0.3, c=seaborn.color_palette("deep", 3)[2], s=35)
plt.plot(realduration, dfit.predict(np.array(realduration).reshape(len(realduration), 1)), linewidth=3)

plt.tight_layout()

# <codecell>

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


plt.figure()
plt.subplot(211)
plt.scatter(s0plotx, s0ploty, c = seaborn.color_palette("deep", 3)[0])
plt.plot(s0plotx, s0.predict(np.array(s0plotx).reshape(len(s0plotx), 1)))
plt.title("S0 vs Real Epidemic Size, Gradient = %f, Intercept = %f" % (s0.coef_[0], s0.predict(0)))
# for j in allss :
   #     plt.scatter(Sbar + Z[e], j[i], c = seaborn.color_palette("deep", 3)[2])
#    print i

plt.subplot(212)
plt.scatter(s0plotx2, s0ploty2, c = seaborn.color_palette("deep", 3)[0], alpha=0.3)
plt.plot(s0plotx2, s1.predict(np.array(s0plotx2).reshape(len(s0plotx2), 1)))
plt.title("S0 vs Simulated Epidemic Size, Gradient = %f, Intercept = %f" % (s1.coef_[0], s1.predict(0)))

plt.tight_layout()
plt.show()

# <codecell>

C.shape

# <codecell>


# <codecell>

plt.

# <codecell>


