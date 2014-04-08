## IMPORTS

# Numerical packages and methods
import numpy as np
from sklearn import linear_model, gaussian_process
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
sensitivity = 5
periodicity = 24
penalty = 1e-3
delay = 8
vaccine = 1965
numSvals = 500
Alpha = 0.97 # set to None to infer
gam = 0. # set to None to infer









# Define epidemics from a time-series
# This function is for the original time-series
def breakepis(x, sensitivity) :

    x2 = x.copy()
    x2[x2 < sensitivity] = 0
    x2 = np.convolve(x2, np.hanning(9), "same")

    z = np.where(x2 > 0)[0]#(np.where(x > sensitivity) and np.where(x2 > sensitivity))[0]  # Find epidemics over sensitivity threshold
    dz = np.where(np.append(np.insert(np.diff(z), 0, 0), -1) != 1)[0]
    
    epi = []
    s = []
    d = []

    for i in range(len(dz) - 1) :
        epi.append(z[dz[i]+3:dz[i+1]-3])
    
    for i, e in enumerate(epi) :
        s.append(np.sum(x[e]))
        d.append(len(e))

    d = np.array(d)
    epi = np.delete(epi, np.where(d == 1)[0])
    s = np.delete(s, np.where(d == 1)[0])
    d = np.delete(d, np.where(d == 1)[0])

    z = []
    for e in epi :
        z = np.hstack((z, e))
            
    return (np.array(epi), s, d, list(z))







# Define epidemics for simulated series
# This function constrains epi durations by the epidemics in the original series
def breaksims(x, sensitivity, realepi) :

    starts = [e[0] for e in realepi]
    starts.append(len(x))

    epi = []
    s = []
    d = []

    for i in range(len(starts)-1) :
        try :
            d.append(np.where(x[starts[i] : starts[i+1]] <= sensitivity)[0][0])
        except IndexError :
            d.append(starts[i+1] - starts[i])

        s.append(np.sum(x[starts[i] : starts[i]+d[-1]]))
        epi.append(range(starts[i], starts[i]+d[-1]))

    return (epi, s, d)







# Bootstrap confidence intervals using time-series sampling
"""

# Bootstrapping
def bootstrap(data, M, statistic) :

    # M times, sample with replacement
    means = np.zeros((M, data.shape[1]))

    for i in range(M) :
        means[i, :] = statistic(data[np.random.randint(0, data.shape[0], 100), :], axis=0)

    stat = np.sort(means, axis=0)

    return (stat[int((alpha/2.0)*M)], stat[int((1-alpha/2.0)*M)])

    #num_samples, statistic, alpha) :
    #n = len(data)
    #idx = np.random.randint(0, n, (num_samples, n))
    #samples = data[idx]
"""






# For ensured full-domain sampling
def pad(X, Y) :
    return np.append(np.insert(X, 0, Y[0]), Y[-1])








# Selective downsampling of cumulative cases
def downsample(X, Y) :
    # Interpolate over even grid
    #X = np.convolve(X, np.hanning(21), "same")
    x = pad(np.linspace(X[0], X[-1], 5*len(X)), X)
    y = pad(interp.interp1d(X, Y)(x), Y)

    # Calculate third derivative
    dy = np.diff(y, n = 3)
    dy[np.abs(dy) < 1] = 0

    # Find zero-crossings
    d  = np.where((dy[:-1] > 0) * (dy[1:] <= 0))[0]
    d = pad(d, range(len(x)))

    return (x[d].reshape(len(x[d]),1), y[d])












def derivative(X, Y) :

    # Central finite difference methods to eighth order coefficients
    c = np.array([1./280, -4./105, 1./5, -4./5, 0, 4./5, -1./5, 4./105, -1./280])

    # Supersample
    x = np.linspace(X[0], X[-1], 10000)
    y = interp.interp1d(X, Y, "linear")(x)
    y[0] = 2*y[1] - y[2]

    # Compute derivative using c
    dx = np.diff(x)[0]
    dy = np.zeros(len(y)-8)

    for i in range(4, len(y)-4) :
        dy[i-4] = np.sum(c * y[i-4:i+5]) / dx


    # Fill in missing left and right values, for now as a stretch
    return interp.interp1d(np.linspace(X[0], X[-1], len(dy)), dy)(X)


























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
    Nt = data["population"].values[:v]


    # Find epidemics
    #z = np.where(C > sensitivity)[0]
    epi, reals, reald, z = breakepis(C, sensitivity)


    # Plot
    plt.figure(figsize=(16, 9), dpi=600)

    plt.subplot(211)
    plt.plot(t, C, linewidth=2)
    plt.title("Reported Cases, %s" % names[idx])
    for e in epi :
        plt.axvspan(t[e[0]], t[e[-1]], color = seaborn.color_palette("deep", 3)[2], alpha=0.3)
    plt.xlabel("Time")
    plt.ylabel("Cases")

    plt.subplot(212)
    plt.plot(t, B, linewidth=2)
    plt.title("Live Births, %s" % names[idx])
    plt.xlabel("Time")
    plt.ylabel("Births")

    plt.tight_layout()
    plt.savefig(prefix + "results/%s_0_timeseries.png" % names[idx])
    print "%s Time-Series done." % names[idx]

    plt.close()














    # Susceptible Reconstruction

    # Compute cumulative births and incidence
    Y = np.cumsum(B)
    X = np.cumsum(C) 

    # Downsample the cumulative plot
    #x, y = downsample(X, Y)
    x = np.linspace(X[0], X[-1], len(X))
    y = interp.interp1d(X, Y)(x)
    y[0] = y[1] - (y[2] - y[1])
    #x = x[:-1].reshape(len(x)-1,1)

    # Rho
    Yhat = gaussian_process.GaussianProcess(nugget = 1e-4)
    Yhat.fit(x.reshape(len(x), 1), y)
    Yhat = Yhat.predict(X.reshape(len(X), 1))
    
    #Yhat = lowess.lowess(y, x.squeeze(), 0.5, return_sorted = False)
    #Yhat = interp.interp1d(x, np.insert(Yhat, 0, Yhat[0]))(X)
    #Yhat = lowess.lowess(Y, X.squeeze(), 0.5, return_sorted = False)


    # Rho : put some splines through Yhat and take the derivative
    rho = derivative(X, Yhat)
    #interp.UnivariateSpline(x, Yhat.predict(x.reshape(len(x),1))).derivative()(np.linspace(X[0], X[-1], len(X)))


    Z = Y - Yhat#.predict(X.reshape(len(X),1))


    # Plots
    plt.figure(figsize=(16, 9), dpi=600)

    plt.subplot(221)
    plt.plot(t, X, linewidth=2)
    plt.plot(t, Y, linewidth=2)
    plt.plot(t, Yhat, linewidth=2)#.predict(X.reshape(len(X), 1)), linewidth=2)
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

    plt.tight_layout()
    plt.savefig(prefix + "results/%s_1_susceptible_reconstruction.png" % names[idx])
    print "%s Susceptible Reconstruction done." % names[idx]
    
    plt.close()






























    if np.any(1./rho > 1) and np.mean(rho) < 1 :
        print "Rho has >1 values; moving onto next geography."
        continue


    # Fit Sbar

    # All possible values of Sbar
    Svals = np.linspace(1, np.abs(np.min(Z))*10, numSvals)

    # Likelihood of fit
    l = np.zeros(len(Svals))



    # Define our parameters
    params = lmfit.Parameters()
    if Alpha is None :
        params.add("alpha", min=0.5, max=.9999, value=0.95) # Alpha

    if gam is None :
        params.add("gamma", min=0.0, max=1.0, value = 0.2)

    for i in range(periodicity) : # Seasonalities
        params.add("r%d" % i, value=0.)

    rstr = ["r%d" % (i % periodicity) for i in z[1:]]

        
    # Objective function
    def profile_residuals(params, rho, C, Z, z, Alpha, gamma, Sestimate) :
        c = C.copy()
        c[np.intersect1d(np.where(C == 0)[0], z).astype(int)] = 1
        
        if Alpha is None : 
            alphafit = params["alpha"].value

        if gamma is None :
            gamma = params["gamma"].value

        r = [params[i].value for i in rstr]
        
        if np.isnan(Sestimate) :
            Sestimate = params["Sest"].value
        
        if Alpha is None :
            return gamma * np.log(Nt[z[:-1]]) + alphafit * np.log(rho[z[:-1]]*c[z[:-1]]) + r + np.log(Sestimate + Z[z[:-1]]) - np.log(rho[z[1:]]*c[z[1:]])

        else :
            return gamma * np.log(Nt[z[:-1]]) + Alpha * np.log(rho[z[:-1]]*c[z[:-1]]) + r + np.log(Sestimate + Z[z[:-1]]) - np.log(rho[z[1:]]*c[z[1:]])
        
        
        


    pbar = progressbar.ProgressBar(widgets = \
            [progressbar.FormatLabel("Evaluating Sbar Likelihoods"), progressbar.Percentage(), progressbar.Bar(), progressbar.ETA()], \
            maxval = numSvals)
    pbar.start()
        
    # Compute best fit for each possible Sbar
    for i, Sestimate in enumerate(Svals) :
        l[i] = lmfit.minimize(profile_residuals, params, args=(rho, C, Z, list(z), Alpha, gam, Sestimate), method="leastsq").chisqr
        pbar.update(i)

    pbar.finish()

    # Purge NaN
    Svals = np.delete(Svals, np.where(np.isnan(l)))
    l = np.delete(l, np.where(np.isnan(l)))
        
    # Fit window
    #fitwindow = 15
    #fitwindowL = np.min([fitwindow, np.argmin(l)])
    #fitwindowR = np.min([fitwindow, len(Svals) - np.argmin(l)])
        
    # Run again using scan estimate
    params.add("Sest", value = Svals[np.argmin(l)])
    L = lmfit.minimize(profile_residuals, params, args=(rho, C, Z, z, Alpha, gam, np.nan), method="leastsq")


    # Extract parameters and errors
    Sbar = L.params["Sest"].value
    r = np.exp([L.params["r" + str(i)].value for i in range(periodicity)])
    alphaSbar = L.params["alpha"].value if Alpha is None else Alpha
    gamma = L.params["gamma"].value if gam is None else gam
    errup = np.exp(np.log(r) + [2*L.params["r" + str(i)].stderr for i in range(periodicity)])
    errdn = np.exp(np.log(r) - [2*L.params["r" + str(i)].stderr for i in range(periodicity)])
        
        
    # Plot
    plt.figure(figsize=(16, 9), dpi=600)

    plt.subplot(121)
    plt.axvline(x=Sbar, color="red", linewidth=2)
    plt.loglog(Svals, l, linewidth=2)
    plt.title("Goodness of Fit, %s" % names[idx])
    plt.xlabel(r"$\bar{S}$")
    plt.ylabel(r"$\chi^2$")
    plt.legend([r"$\bar{S}$ = %d, $\alpha$ = %.03f, $\gamma$ = %.03f" % (Sbar, alphaSbar, gamma)])


    plt.subplot(122)
    plt.plot(r, linewidth=3)
    plt.fill_between(range(periodicity), errup, errdn, color=colours[0], alpha=0.3)
    plt.xlim([0, periodicity])
    plt.title("Periodicity")
    plt.xlabel("Period")
    
    plt.tight_layout()
    plt.savefig(prefix + "results/%s_2_meanS_periodicity.png" % names[idx])
    print "%s Mean S and Periodicity done." % names[idx]

    plt.close()


























    # Simulations

    allss = []
    allsd = []
    allrs = []
    allrd = []
    allepi = []
    alliei = []
    I = []
    pbar = progressbar.ProgressBar(widgets = \
            [progressbar.FormatLabel("Running Simulations"), progressbar.Percentage(), progressbar.Bar(), progressbar.ETA()], \
            maxval = int(sys.argv[2]))
    pbar.start()

    for q in range(int(sys.argv[2])) :
        
        predI = np.zeros_like(C)
        predS = np.zeros_like(C)
        starts = [e[0] for e in epi]

        
        # Seed initial epidemic points
        for index, e in enumerate(starts) :
            ss = []
            predI[e] = np.ceil(rho[e] * C[e])
            predS[e] = np.ceil(Sbar + Z[e])
        
            for i in range(e+1, len(C)) :
                
                """
                predI[i] = np.random.negative_binomial(max(np.round(predI[i-1])**alphaSbar, 1), \
                                                       max(np.round(predI[i-1])**alphaSbar, 1) / ( max(np.round(predI[i-1]), 1) + \
                                                                   r[i % periodicity] * ( predI[i-1] ** alphaSbar ) * predS[i-1] / (Nt[i-1]**gamma)))
                """
                bsi = r[i % periodicity] * predS[i-1] * (predI[i-1] ** (alphaSbar-1)) if np.isfinite(predI[i-1] ** (alphaSbar-1)) else 0

                predI[i] = np.random.negative_binomial(max(np.round(predI[i-1]), 1), \
                                                       1. / ( 1. + bsi ))

                predS[i] = max(B[max(i - delay, 0)] + predS[i-1] - predI[i], 0)

            
            #simsizes.append(np.sum(predI[e:]))
            #simduration.append(min(np.argmin(predI[e:]), maxdur[index]))
            #realduration.append(np.argmin(C[e:] * rho[e:]))
            #realsizes.append(np.sum(rho[e:e+realduration[-1]] * C[e:e+realduration[-1]]))



        sepi, simsize, simdur = breaksims(predI, sensitivity, epi)

        allss.append(simsize)
        allrs.append(reals)
        allsd.append(simdur)
        allrd.append(reald)
        allepi.append(sepi)
        alliei.append(np.diff(starts))# - simdur[:-1])
        I.append(predI)
        pbar.update(q)

    pbar.finish()

        
    allieisize = np.array([s[1:] for s in allss]).ravel()
    allieidur = np.array([d[1:] for d in allsd]).ravel()
    allss = np.array(allss).ravel()
    allrs = np.array(allrs).ravel()
    allsd = np.array(allsd).ravel()
    allrd = np.array(allrd).ravel()
    allepi = np.array(allepi).ravel()
    alliei = np.array(alliei).ravel()
        
    #idx = allrs > 500
    #allss = allss[idx].reshape(np.sum(idx), 1)
    #allrs = allrs[idx].reshape(np.sum(idx), 1)

     
    sslope, sintercept, sr, sp, _ = st.linregress(allrs.squeeze(), allss.squeeze())
    dslope, dintercept, dr, dp, _ = st.linregress(allrd.squeeze(), allsd.squeeze())

    xs, xd = np.linspace(0, allrs.max(), 500), np.linspace(0, allrd.max(), 500)
    ys, yd = sslope * xs + sintercept, dslope * xd + dintercept
    #sfit = linear_model.BayesianRidge(fit_intercept=False)
    #dfit = linear_model.BayesianRidge(fit_intercept=False)

    #dfit.fit(allrd, allsd)
    #sfit.fit(allrs, allss)

    I = np.array(I)    
        


    # Boostrap confidence intervals
    low = np.zeros(I.shape[1])
    high = np.zeros(I.shape[1])

    """pbar = progressbar.ProgressBar(widgets = \
            [progressbar.Percentage(), progressbar.Bar(), progressbar.ETA()], \
            maxval = int(I.shape[1]))
    #pbar.start()"""

#   for i in range(I.shape[1]) :
#       low[i], high[i] = bootstrap(I[:, i], 1000, np.mean, 0.95)
#       pbar.update(i)

        
    # Plot   
    
    plt.figure(figsize=(16, 9), dpi=600)
#   plt.fill_between(t, low, high, color = colours[2], linewidth=1, alpha=0.4)
    plt.plot(t, np.mean(I, axis=0), color = colours[2], linewidth=2)
    plt.plot(t, C*rho, c = colours[0], linewidth=2, alpha = 0.8)

    plt.tight_layout()
    plt.savefig(prefix + "results/%s_3_predictions.png" % names[idx])
    print "%s Predictions done." % names[idx]

    plt.close()



    # Size and Duration of Epidemics : Real vs Predicted
    plt.figure(figsize=(16, 9), dpi=600)

    plt.subplot(211)
    plt.title("%s, Sizes : Slope = %.3f, Intercept = %.1f, R^2 = %.3f, p = %e" % (names[idx], sslope, sintercept, sr, sp))
    plt.xlabel("Real Size")
    plt.ylabel("Simulated Size")
    nepi = len(allrs) / (q+1)
    plt.errorbar(allrs[:nepi], np.mean(allss.reshape(q+1, nepi), axis=0), 2*np.std(allss.reshape(q+1, nepi), axis=0), fmt="o", ms=15, c=colours[2])
    for i in range(q) :
        plt.scatter(allrs[i], allss[i], alpha=0.3, c=colours[2], s=35)
    plt.plot(xs, ys, linewidth=2, c=colours[0])
        
    plt.subplot(212)
    plt.title("Durations : Slope = %.3f, Intercept = %.1f, R^2 = %.3f, p = %e" % (dslope, dintercept, dr, dp))
    plt.xlabel("Real Duration")
    plt.ylabel("Simulated Duration")
    plt.errorbar(allrd[:nepi], np.mean(allsd.reshape(q+1, nepi), axis=0), 2*np.std(allsd.reshape(q+1, nepi), axis=0), fmt="o", ms=15, c=colours[2])
    for i in range(q) :
        plt.scatter(allrd[i], allsd[i], alpha=0.3, c=colours[2], s=35)
    plt.plot(xd, yd, linewidth=2)

    plt.tight_layout()
    plt.savefig(prefix + "results/%s_4_sizes_durations.png" % names[idx])
    print "%s Sizes and Durations done." % names[idx]

    plt.close()



















    # Susceptibles vs Sizes 

    allss = np.array(allss).ravel()

    s0 = np.array([Sbar + Z[e[0]] for e in epi])#np.array([np.mean(Sbar + Z[e]) for e in epi])
    slopes0, intercepts0, rs0, ps0, _ = st.linregress(s0[reals > 20], reals[reals > 20])

    s1 = np.array([np.mean(Sbar + Z[e]) for e in np.array(allepi).ravel()])
    slopes1, intercepts1, rs1, ps1, _ = st.linregress(s1[allss > 20], allss[allss > 20])

    s0x = np.linspace(0, s0.max(), 500)
    s0y = s0x * slopes0 + intercepts0

    s1x = np.linspace(0, s1.max(), 500)
    s1y = s1x * slopes1 + intercepts1

    

    # Plot
    plt.figure(figsize=(16, 9), dpi=600)

    plt.subplot(211)
    plt.scatter(s0, reals, c = colours[0])
    plt.plot(s0x, s0y, linewidth=2)
    plt.title("%s S0 vs Real Size, Slope = %.3f, Intercept = %.1f, R^2 = %.3f, p = %e" % (names[idx], slopes0, intercepts0, rs0, ps0))


    plt.subplot(212)
    plt.scatter(s1, allss, c = colours[0], alpha=0.3)
    plt.plot(s1x, s1y, linewidth=2)
    plt.title("S0 vs Simulated Size, Slope = %.3f, Intercept = %.1f, R^2 = %.3f, p = %e" % (slopes1, intercepts1, rs1, ps1))

    plt.tight_layout()
    plt.savefig(prefix + "results/%s_5_s0_vs_size.png" % names[idx])
    print "%s S0 vs Sizes done." % names[idx]

    plt.close()






















    # Inter-epidemic intervals
    plt.figure(figsize=(16, 9), dpi=600)

    ieix = np.linspace(alliei.min(), alliei.max(), 500)
    rieiss, rieisi, rieisr, rieisp, _ = st.linregress(np.diff(starts), reals[1:])
    rieids, rieidi, rieidr, rieidp, _ = st.linregress(np.diff(starts), reald[1:])
    ieiss, ieisi, ieisr, ieisp, _ = st.linregress(alliei, allieisize)
    ieids, ieidi, ieidr, ieidp, _ = st.linregress(alliei, allieidur)

    plt.subplot(211)
    plt.scatter(alliei, allieisize, alpha=0.2, c=colours[0])
    plt.scatter(np.diff(starts), reals[1:], s=100, c=colours[2])
    plt.plot(ieix, rieiss * ieix + rieisi, linewidth=2)
    plt.plot(ieix, ieiss * ieix + ieisi, linewidth=2)
    plt.title("%s, IEI vs Size : Slope = %.3f, Intercept = %.1f, R^2 = %.3f, p = %e" % (names[idx], ieiss, ieisi, ieisr, ieisp))
    plt.xlabel("Interepidemic Interval (biweeks)")
    plt.ylabel("Size of Epidemic")
    plt.legend(["Real Fit", "Sim Fit", "Simulated", "Real"])

    plt.subplot(212)
    plt.scatter(alliei, allieidur, alpha=0.2, c=colours[0])
    plt.scatter(np.diff(starts), reald[1:], s=100, c=colours[2])
    plt.plot(ieix, rieids * ieix + rieidi, linewidth=2)
    plt.plot(ieix, ieids * ieix + ieidi, linewidth=2)
    plt.title("IEI vs Duration : Slope = %.3f, Intercept = %.1f, R^2 = %.3f, p = %e" % (ieids, ieidi, ieidr, ieidp))
    plt.xlabel("Interepidemic Interval (biweeks)")
    plt.ylabel("Duration of Epidemic")
    plt.legend(["Real Fit", "Sim Fit", "Simulated", "Real"])

    plt.tight_layout()
    plt.savefig(prefix + "results/%s_6_iei.png" % names[idx])
    print "IEI done."














    # R estimates
    plt.figure(figsize=(16,9), dpi=600)

    # Starting points of each epidemic
    L = []
    for e in epi[:-1] :
        L.append(len(e))
    L = np.cumsum(L)

    # R vs t
    Reff = C[1:].astype(float) / C[:-1]

    # S0 at beginning of each epi, and max R for each of those epis
    S0 = [Sbar + Z[e[0]] for e in epi]
    Rm = [Reff[e].max() for e in epi]



    plt.subplot2grid((2,2), (0,0))
    plt.plot(C[z], linewidth=2)
    for e in L :
        plt.axvline(e, c=colours[2])
    plt.title("%s, Incidence" % names[idx])
    plt.ylabel("Cases")

    plt.subplot2grid((2, 2), (1, 0))
    plt.plot(Reff[z], linewidth=2)
    plt.axhline(1, c=colours[1])
    for e in L :
        plt.axvline(e, c=colours[2])
    plt.title("Effection Reproduction Ratio")
    plt.ylabel("$R_{eff}$")

    plt.subplot2grid((2,2), (0, 1), rowspan=2)
    plt.scatter(S0, Rm)
    plt.title("$R_0$ as Max $R_{eff}$ vs $S_0$")
    plt.xlabel("$S_0$")
    plt.ylabel("max[$R_eff$]")

    plt.tight_layout()
    plt.savefig(prefix + "results/%s_7_r0.png" % names[idx])
    print "R0 done."


