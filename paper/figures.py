 # -*- coding: latin-1

# Generate figures for paper

# Plots to generate :
#	Time-series / Predictions - alpha
#	Reporting rates
#	Periodicity
#	Predictions
#	Size

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import matplotlib as mpl
import seaborn
#mpl.rcParams['text.latex.preamble'] = [r'\usepackage[utf8]{fontenc}'] 
mpl.rcParams['text.latex.unicode'] = True 
mpl.rcParams['text.usetex'] = True   

colours = seaborn.color_palette("deep", 8)
scalefactor = 20.
xdim = 1.1
ydim = 1.


name = ["Akureyri", "Bornholm", "Faroe Islands", r"Hafnarfj\"{o}rdur", r"Reykjav\'{i}k", "Vestmannaeyjar"]
#name = []
t = []
ts = []
pred = []
rho = []
r = []
rup = []
rdn = []
sizex = []
sizey = []
sizeerrx = []
sizeerry = []
sizeerre = []
r2 = []
p = []
pearson = []
pearsonzero = []
sbar = []
sn = []
Z = []
grad = []
ciu = []
cid = []
alpha = []
b = []






for file in os.listdir("figures/") :
	if file.endswith(".json") :
		data = pd.read_json("figures/%s" % file)
		#name.append(" %s" % file.split(".json")[0])
		t.append(data["t"].values[0])
		ts.append(data["ts"].values[0])
		pred.append(data["pred"].values[0])
		rho.append(data["rho"].values[0])
		r.append(np.array(data["r"].values[0]))
		rup.append(np.array(data["rup"].values[0]))
		rdn.append(np.array(data["rdn"].values[0]))
		sizex.append(data["sizex"].values[0])
		sizey.append(data["sizey"].values[0])
		sizeerrx.append(data["sizeerrx"].values[0])
		sizeerry.append(data["sizeerry"].values[0])
		sizeerre.append(data["sizeerre"].values[0])
		r2.append(data["r2"].values[0])
		p.append(data["p"].values[0])
		pearson.append(data["pearson"].values[0])
		pearsonzero.append(data["pearsonzero"].values[0])
		sbar.append(data["sbar"].values[0])
		grad.append(data["grad"].values[0])
		sn.append(data["sn"].values[0])
		ciu.append(data["predciu"].values[0])
		cid.append(data["predcid"].values[0])
		alpha.append(data["alpha"].values[0])
		b.append(data["b"].values[0])




print "ALPHA : ", alpha



# Births
B = []
for file in os.listdir("../data/iceland") :
	if file.endswith(".csv") :
		if file != "all.csv" :
			B.append(pd.read_csv("../data/iceland/" + file)["births"].values)

B.insert(1, pd.read_csv("../data/bornholm/bornholm.csv")["births"].values)
B.insert(2, pd.read_csv("../data/faroe/faroe.csv")["births"].values)



for i in range(len(B)) :
	v = np.where(np.array(t[i]) > 1965)[0][0] if t[i][-1] > 1965 else len(t[i])
	B[i] = B[i][:v]





fig, ax = plt.subplots(len(name), 1, sharex=True, figsize=(xdim*210/scalefactor, ydim*297/scalefactor))
for i in range(len(name)) :

	ax[i].plot(t[i], ts[i], lw=2, c=colours[0])
	ax[i].set_title("%s" % name[i].decode("utf-8"), fontsize=16, loc="left")
	ax[i].tick_params(labelsize=16)
	ax[i].locator_params(nbins=5, axis="y")
	ax[i].set_xlim([1901, 1965])
	ax[i].set_ylim([0, np.max(ts[i])*1.1])

	#ax[i][1].plot(t[i], B[i], lw=2, c=colours[2])
	#ax[i][1].tick_params(labelsize=16)
	#ax[i][1].locator_params(nbins=5, axis="y")
	#ax[i][1].set_xlim([np.min(t[i]), np.max(t[i])])
	
	#ax2 = ax1.twinx()
	#ax2.plot(t[i], rho[i], lw=3, c=colours[2], alpha=0.6)
	#ax2.grid(False)
	#ax2.tick_params(labelsize=16)
	#ax2.locator_params(nbins=5, axis="y")
	

	
	#ax2.set_xlim([np.min(t[i]), np.max(t[i])])

	#if i == 0 :
	#	ax[i].legend(["Incidence", "Reporting Rate"], loc=2)

plt.figtext(.015, 0.5, r"Reported Incidence, $C_t$", ha="center", va="center", rotation="vertical", fontsize=16)
#plt.figtext(.99, 0.5, r"Reporting Rate, $1/\rho$", ha="center", va="center", rotation="vertical", fontsize=16)
plt.figtext(.55, 0.01, "Time", ha="center", va="center", fontsize=16)
plt.tight_layout(rect=(0.015, 0.015, .98, 1))
plt.savefig("figures/0_incidence.pdf")
plt.close()







"""
fig, axs = plt.subplots(len(name), 2, sharex=True, figsize=(xdim*210/scalefactor, ydim*297/scalefactor))
##plt.suptitle("Seasonality", fontsize=20, y=.999)
for i in range(len(name)) :
	axs[i][0].plot(r[i] * sbar[i], lw=2)
	axs[i][0].fill_between(range(24), rup[i] * sbar[i], rdn[i] * sbar[i], color=colours[0], alpha=0.3)
	axs[i][0].set_title(name[i], fontsize=16, loc="left")
	axs[i][0].set_xlim([0, 23])
	axs[i][0].set_xticks(range(0, 23, 4))
	axs[i][0].tick_params(labelsize=16)
	axs[i][0].locator_params(nbins=5, axis="y")

	axs[i][1].plot(t[i], rho[i], lw=2)
	axs[i][1].grid(False)
	axs[i][1].tick_params(labelsize=16)
	axs[i][1].locator_params(nbins=5, axis="y")
	axs[i][1].set_xlim([np.min(t[i]), np.max(t[i])])


plt.figtext(.5, 0.01, "Period", ha="center", va="center", fontsize=16)
plt.figtext(.01, 0.5, r"Seasonality, $r \bar{S}$", ha="center", va="center", rotation="vertical", fontsize=16)
plt.tight_layout(rect=(0.015, 0.015, .99, 1)) # (left, bottom, right, top) 
plt.savefig("figures/1_seasonality.pdf")
plt.close()
"""





"""
fig, axs = plt.subplots(len(name), 2, sharex='col', figsize=(xdim*210/scalefactor, ydim*297/scalefactor))
##plt.suptitle("Predicted Cases", fontsize=20, y=.999)
for i in range(len(name)) :
	axs[i][0].plot(t[i], ts[i], lw=2, c=colours[0])
	axs[i][0].plot(t[i], pred[i], lw=1, c=colours[2], alpha=1)
	axs[i][0].fill_between(t[i], cid[i], ciu[i], color=colours[2], alpha=0.3)
	axs[i][0].tick_params(labelsize=16)
	axs[i][0].locator_params(nbins=5, axis="y")
	axs[i][0].set_title(u"%s. Zero-corrected $R^2$ = %.03f" % (name[i].decode("utf-8"), pearsonzero[i]), fontsize=16, loc="left")
	axs[i][0].set_xlim([1901, 1965])
	axs[i][0].set_ylim([0, np.max([np.max(ts[i]), np.max(pred[i])])*1.1])
	if i == 0 :
		axs[i][0].legend([r"Observed Incidence", "Predicted Incidence"], loc=2)


	axs[i][1].plot(r[i] * sbar[i], lw=2)
	axs[i][1].fill_between(range(24), rup[i] * sbar[i], rdn[i] * sbar[i], color=colours[0], alpha=0.3)
	#axs[i][1].set_title(name[i], fontsize=16, loc="left")
	axs[i][1].set_xlim([0, 23])
	axs[i][1].set_xticks(range(0, 23, 4))
	axs[i][1].tick_params(labelsize=16)
	axs[i][1].locator_params(nbins=5, axis="y")

plt.figtext(.015, 0.5, r"Observed Incidence, $C_t$ and Predicted Incidence, $I_t / \rho_t$", ha="center", va="center", rotation="vertical", fontsize=16)
plt.figtext(.99, 0.5, r"Seasonality, $r_t\,\bar{S}$", ha="center", va="center", rotation="vertical", fontsize=16)
plt.figtext(.285, 0.011, "Time", ha="center", va="center", fontsize=16)
plt.figtext(.751, 0.011, "Period", ha="center", va="center", fontsize=16)
plt.tight_layout(rect=(0.02, 0.02, .97, 1))
plt.savefig("figures/2_predictions_.pdf")
plt.close()
"""






"""
fig, axs = plt.subplots(len(name), 1, figsize=(xdim*210/scalefactor, ydim*297/scalefactor))
#fig.suptitle("Predicted Epidemic Sizes", fontsize=20, y=.999)
for i in range(len(name)) :
	axs[i].errorbar(sizeerrx[i], np.array(sizeerry[i]) / np.array(sizeerrx[i]), yerr = np.array(sizeerre[i]) / np.array(sizeerrx[i]), fmt="o", ms=8, c=colours[0])
	#axs[i].set_yscale("log")
	#axs[i].plot(sizex[i], np.array(sizey[i]) / np.array(sizex[i]), lw=2, c=colours[2])
	axs[i].axhline(np.mean(np.array(sizeerry[i]) / np.array(sizeerrx[i])), color=colours[2], lw=2)
	axs[i].tick_params(labelsize=16)
	axs[i].locator_params(nbins=5, axis="y")
	axs[i].set_xlim([0, np.max(sizex[i])*1.02])
	axs[i].set_title(r"%s. $\bar{q}$ = %.02f" % (name[i].decode("utf-8"), grad[i]), fontsize=16, loc="left")

plt.figtext(.01, 0.5, "Predicted Size Coefficient $q$", ha="center", va="center", rotation="vertical", fontsize=16)
plt.figtext(.5, 0.01, "Actual Epidemic Size", ha="center", va="center", fontsize=16)
plt.tight_layout(rect=(0.01, 0.015, .99, 1))
plt.savefig("figures/3_sizes.pdf")
plt.close()
"""







fig, axs = plt.subplots(3,2, figsize=(xdim*210/scalefactor, ydim*180/scalefactor))
#fig.suptitle("Predicted Epidemic Sizes", fontsize=20, y=.999)
for i in range(len(name)) :
	axs[i%3][int(i>2)].plot(sizex[i], np.array(sizey[i]), lw=2, c=colours[2])
	axs[i%3][int(i>2)].plot([0, np.max(sizex[i])], [0, np.max(sizex[i])], lw=2, c=colours[1], alpha=0.7)
	axs[i%3][int(i>2)].errorbar(sizeerrx[i], np.array(sizeerry[i]), yerr = np.array(sizeerre[i]), fmt="o", ms=8, c=colours[0])
	#plt.axhline(1, color=colours[2], lw=2)
	axs[i%3][int(i>2)].tick_params(labelsize=16)
	axs[i%3][int(i>2)].set_ylim([0, 1.15*np.max([np.max(sizeerry[i]), np.max(sizey[i])])])
	axs[i%3][int(i>2)].locator_params(nbins=5, axis="y")
	axs[i%3][int(i>2)].set_xlim([0, np.max(sizex[i])*1.02])
	axs[i%3][int(i>2)].text(0.04, 0.85, "%s." % name[i].decode("utf-8"), fontsize=16, horizontalalignment="left", transform=axs[i%3][int(i>2)].transAxes)
	#axs[i%3][int(i>2)].text(0.04, 0.85, "%s. $R^2$ = %.2f, gradient = %.02f" % (name[i].decode("utf-8"), r2[i], grad[i]), fontsize=16, horizontalalignment="left", transform=axs[i%3][int(i>2)].transAxes)

	#if i == 0 :
	#	axs[i%3][int(i>2)].legend(["Regression", "Identity"], loc=2, prop={"size" : 14})

plt.figtext(.015, 0.5, r"Predicted Epidemic Size, $I_t / \rho_t$", ha="center", va="center", rotation="vertical", fontsize=16)
plt.figtext(.5, 0.01, r"Observed Epidemic Size, $C_t$", ha="center", va="center", fontsize=16)
plt.tight_layout(rect=(0.03, 0.02, .99, 1))
plt.savefig("figures/q3.pdf")
plt.close()
print "Figure 3 made."





"""
plt.figure(figsize=(xdim*210/scalefactor, ydim*297/scalefactor))
for i in range(len(name)) :

	T = np.array(ts[i]) if len(ts[i]) % 2 == 0 else np.array(ts[i][:-1])
	T = T.reshape(len(T)/2, 2).sum(axis=1)

	num = np.zeros(12)
	count = np.zeros(12)

	for j, TS in enumerate(ts[i]) :
		if TS > 0.5 :
			num[j % 12] += 1
		count[j % 12] += 1

	plt.subplot(len(name), 1, i+1)
	plt.plot(num / count.astype(float), lw=2)
	plt.tick_params(labelsize=16)
	plt.locator_params(nbins=5, axis="y")

plt.figtext(.01, 0.5, "Fraction of Month with Cases", ha="center", va="center", rotation="vertical", fontsize=16)
plt.figtext(.5, 0.01, "Month", ha="center", va="center", fontsize=16)
plt.tight_layout(rect=(0.02, 0.015, .99, 1))
plt.savefig("figures/5_infectiontiming.pdf")
plt.close()
"""













fig, axs = plt.subplots(3, 2, sharex=True, figsize=(xdim*210/scalefactor, ydim*115/scalefactor))
##plt.suptitle("Seasonality", fontsize=20, y=.999)
for i in range(len(name)) :
	axs[i % 3][np.floor(i/3)].plot(t[i], ts[i], lw=2)
	axs[i % 3][np.floor(i/3)].set_title("%s" % name[i].decode("utf-8"), fontsize=16, loc="left")
	axs[i % 3][np.floor(i/3)].set_xlim([1900, 1965])
	#axs[i % 3][np.floor(i/3)].set_xticks(range(0, 23, 4))
	axs[i % 3][np.floor(i/3)].tick_params(labelsize=14)
	axs[i % 3][np.floor(i/3)].locator_params(nbins=5, axis="y")

	a = axs[i%3][np.floor(i/3)].twinx()
	a.plot(t[i], rho[i], lw=3, c=colours[2], alpha=0.6)
	a.grid(False)
	a.tick_params(labelsize=14)
	a.locator_params(nbins=5, axis="y")
	a.set_xlim([1900, 1965])

plt.figtext(.015, 0.5, r"Reported Incidence", ha="center", va="center", rotation="vertical", fontsize=16)
plt.figtext(.5, 0.05, "Time", ha="center", va="center", fontsize=16)
plt.tight_layout(rect=(0.02, 0.05, .99, 1)) # (left, bottom, right, top) 
plt.savefig("figures/6_fortalk.pdf")
plt.close()








fig, axs = plt.subplots(3, 2, sharex=True, figsize=(xdim*210/scalefactor, ydim*115/scalefactor))
##plt.suptitle("Seasonality", fontsize=20, y=.999)
for i in range(len(name)) :
	axs[i % 3][np.floor(i/3)].plot(r[i] * sbar[i], lw=2)
	axs[i % 3][np.floor(i/3)].fill_between(range(24), rup[i] * sbar[i], rdn[i] * sbar[i], color=colours[0], alpha=0.3)
	axs[i % 3][np.floor(i/3)].set_title(name[i].decode("utf-8"), fontsize=16, loc="left")
	axs[i % 3][np.floor(i/3)].set_xlim([0, 23])
	axs[i % 3][np.floor(i/3)].set_xticks(range(0, 23, 4))
	axs[i % 3][np.floor(i/3)].tick_params(labelsize=16)
	axs[i % 3][np.floor(i/3)].locator_params(nbins=5, axis="y")
plt.figtext(.5, 0.05, "Period", ha="center", va="center", fontsize=16)
plt.figtext(.01, 0.5, r"Seasonality, $r \bar{S}$", ha="center", va="center", rotation="vertical", fontsize=16)
plt.tight_layout(rect=(0.02, 0.05, .99, 1)) # (left, bottom, right, top) 
plt.savefig("figures/7_fortalk.pdf")
plt.close()










fig, axs = plt.subplots(3, 2, sharex=True, figsize=(xdim*210/scalefactor, ydim*115/scalefactor))
##plt.suptitle("Predicted Cases", fontsize=20, y=.999)
for i in range(len(name)) :
	axs[i % 3][np.floor(i/3)].plot(t[i], ts[i], lw=2, c=colours[0])
	axs[i % 3][np.floor(i/3)].plot(t[i], pred[i], lw=1, c=colours[2], alpha=1)
	axs[i % 3][np.floor(i/3)].fill_between(t[i], cid[i], ciu[i], color=colours[2], alpha=0.3)
	axs[i % 3][np.floor(i/3)].tick_params(labelsize=16)
	axs[i % 3][np.floor(i/3)].locator_params(nbins=5, axis="y")
	axs[i % 3][np.floor(i/3)].set_title("%s. $R^2$ = %.02f" % (name[i].decode("utf-8"), pearson[i]), fontsize=16, loc="left")
	axs[i % 3][np.floor(i/3)].set_xlim([1900, 1965])
	axs[i % 3][np.floor(i/3)].set_ylim([0, max(np.max(ts[i]), np.max(pred[i]))])
	if i == 0 and np.floor(i/3) == 0 :
		axs[i % 3][np.floor(i/3)].legend([r"Observed Incidence", "Predicted Incidence"], loc=2)

plt.figtext(.015, 0.5, r"Incidence, $\rho_t\,C_t$", ha="center", va="center", rotation="vertical", fontsize=16)
plt.figtext(.5, 0.02, "Time", ha="center", va="center", fontsize=16)
plt.tight_layout(rect=(0.015, 0.05, .99, 1))
plt.savefig("figures/8_fortalk.pdf")
plt.close()











fig, axs = plt.subplots(3, 2, figsize=(xdim*210/scalefactor, ydim*115/scalefactor))
#fig.suptitle("Predicted Epidemic Sizes", fontsize=20, y=.999)
for i in range(len(name)) :
	axs[i % 3][np.floor(i/3)].errorbar(sizeerrx[i], sizeerry[i], yerr = sizeerre[i], fmt="o", ms=8, c=colours[0])
	axs[i % 3][np.floor(i/3)].plot(sizex[i], sizey[i], lw=2, c=colours[2])
	axs[i % 3][np.floor(i/3)].tick_params(labelsize=16)
	axs[i % 3][np.floor(i/3)].locator_params(nbins=5, axis="y")
	axs[i % 3][np.floor(i/3)].set_xlim([0, np.max(sizex[i])*1.02])
	axs[i % 3][np.floor(i/3)].set_ylim([0, max(np.max(sizey[i]), np.max(sizeerry[i]))*1.1])
	axs[i % 3][np.floor(i/3)].set_title("%s. $R^2$ = %.2f, slope = %.02f" % (name[i].decode("utf-8"), r2[i], grad[i]), fontsize=16, loc="left")

plt.figtext(.015, 0.5, "Simulated Epidemic Size", ha="center", va="center", rotation="vertical", fontsize=16)
plt.figtext(.5, 0.02, "Actual Epidemic Size", ha="center", va="center", fontsize=16)
plt.tight_layout(rect=(0.02, 0.035, .99, 1))
plt.savefig("figures/9_fortalk.pdf")
plt.close()


















fig, axs = plt.subplots(3, 2, figsize=(xdim*210/scalefactor, ydim*115/scalefactor))
#fig.suptitle("Predicted Epidemic Sizes", fontsize=20, y=.999)
for i in range(len(name)) :
	axs[i % 3][np.floor(i/3)].scatter(b[i], sizeerrx[i][1:], c=colours[0])
	axs[i % 3][np.floor(i/3)].plot(b[i], b[i], lw=2)
plt.tight_layout(rect=(0.02, 0.035, .99, 1))
plt.savefig("figures/10.pdf")
plt.close()

























fig, axs = plt.subplots(6, 1, sharex=True, figsize=(xdim*210/scalefactor, ydim*180/scalefactor))
##plt.suptitle("Predicted Cases", fontsize=20, y=.999)
for i in range(len(name)) :
	axs[i].plot(t[i], ts[i], lw=2, c=colours[0])
	axs[i].plot(t[i], -np.array(pred[i]), lw=2, c=colours[2], alpha=1)
	axs[i].plot(t[i], ts[i], lw=2, c=colours[0], alpha=0.5)
	axs[i].fill_between(t[i], -np.array(cid[i]), -np.array(ciu[i]), color=colours[2], alpha=0.4)
	axs[i].tick_params(labelsize=16)
	axs[i].locator_params(nbins=4, axis="y")
	axs[i].text(0.019, .2, "%s." % name[i].decode("utf-8"), fontsize=16, horizontalalignment="left", transform=axs[i].transAxes)
	#axs[i].text(0.019, .2, "%s. $R^2$ = %.02f" % (name[i].decode("utf-8"), pearson[i]), fontsize=16, horizontalalignment="left", transform=axs[i].transAxes)
	axs[i].set_xlim([1901, 1965])
	axs[i].set_ylim([-1.1*np.max(pred[i]), 1.1*np.max(ts[i])])
	#if i == 0 :
	#	axs[i].legend([r"Observed Incidence", "Predicted Incidence"], loc=3, prop={"size" : 14})

plt.figtext(.02, 0.5, r"Observed Incidence, $C_t$ and Predicted Incidence, $I_t / \rho_t$", ha="center", va="center", rotation="vertical", fontsize=16)
plt.figtext(.5, 0.008, "Time", ha="center", va="center", fontsize=16)
plt.tight_layout(rect=(0.025, 0.02, .99, 1))
plt.savefig("figures/q1.pdf")
plt.close()
print "Figure 1 made."
















fig, axs = plt.subplots(len(name), 2, sharex='col', figsize=(xdim*210/scalefactor, ydim*180/scalefactor))
##plt.suptitle("Predicted Cases", fontsize=20, y=.999)
for i in range(len(name)) :
	axs[i][0].plot(t[i], rho[i], lw=2)
	axs[i][0].tick_params(labelsize=16)
	axs[i][0].locator_params(nbins=5, axis="y")
	axs[i][0].set_xlim([1901, 1965])
	axs[i][0].locator_params(nbins=5, axis="y")
	axs[i][0].text(0.035, .7, name[i].decode("utf-8"), fontsize=16, horizontalalignment="left", transform=axs[i][0].transAxes)
	axs[i][0].set_ylim([0.85*np.min(rho[i]), 1.15*np.max(rho[i])])


	axs[i][1].plot(r[i] * sbar[i], lw=2)
	axs[i][1].fill_between(range(24), rup[i] * sbar[i], rdn[i] * sbar[i], color=colours[0], alpha=0.3)
	axs[i][1].set_xlim([0, 23])
	axs[i][1].set_xticks(range(0, 23, 4))
	axs[i][1].tick_params(labelsize=16)
	axs[i][1].locator_params(nbins=4, axis="y")

plt.figtext(.015, 0.5, r"Reporting rate $1/\rho_t$", ha="center", va="center", rotation="vertical", fontsize=16)
plt.figtext(.99, 0.5, r"Seasonality, $r_t\,\bar{S}$", ha="center", va="center", rotation="vertical", fontsize=16)
plt.figtext(.285, 0.011, "Time", ha="center", va="center", fontsize=16)
plt.figtext(.751, 0.011, "Period", ha="center", va="center", fontsize=16)
plt.tight_layout(rect=(0.03, 0.025, .97, 1))
plt.savefig("figures/q2.pdf")
plt.close()
print "Figure 2 made."





