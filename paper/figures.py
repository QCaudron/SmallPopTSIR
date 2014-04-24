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
from matplotlib import rc
import seaborn
rc("text", usetex=True)

colours = seaborn.color_palette("deep", 8)
scalefactor = 20.



name = []
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







for file in os.listdir("figures/") :
	if file.endswith(".json") :
		data = pd.read_json("figures/%s" % file)
		name.append(file.split(".json")[0])
		t.append(data["t"].values[0])
		ts.append(data["ts"].values[0])
		pred.append(data["pred"].values[0])
		rho.append(data["rho"].values[0])
		r.append(data["r"].values[0])
		rup.append(data["rup"].values[0])
		rdn.append(data["rdn"].values[0])
		sizex.append(data["sizex"].values[0])
		sizey.append(data["sizey"].values[0])
		sizeerrx.append(data["sizeerrx"].values[0])
		sizeerry.append(data["sizeerry"].values[0])
		sizeerre.append(data["sizeerre"].values[0])
		r2.append(data["r2"].values[0])
		p.append(data["p"].values[0])
		pearson.append(data["pearson"].values[0])







plt.figure(figsize=(210/scalefactor, 297/scalefactor))
plt.suptitle("Reported Incidence and Inferred Reporting Rate", fontsize=14, y=.999)
for i in range(len(name)) :

	ax1 = plt.subplot(len(name), 1, i+1)
	ax1.plot(t[i], ts[i], lw=2)
	ax1.plot(t[i], ts[i], lw=2, c=colours[2], alpha=0.6)
	ax1.plot(t[i], ts[i], lw=2, c=colours[0])
	ax1.set_title(name[i], fontsize=12, loc="left")
	
	ax2 = ax1.twinx()
	ax2.plot(t[i], rho[i], lw=3, c=colours[2], alpha=0.6)
	ax2.grid(False)

	ax1.set_xlim([np.min(t[i]), np.max(t[i])])
	ax2.set_xlim([np.min(t[i]), np.max(t[i])])

	if i == 0 :
		ax1.legend(["Incidence", "Reporting Rate"], loc=2)

plt.figtext(.01, 0.5, r"Reported Incidence, $C$", ha="center", va="center", rotation="vertical")
plt.figtext(.99, 0.5, r"Reporting Rate, $1/\rho$", ha="center", va="center", rotation="vertical")
plt.figtext(.5, 0.01, "Time", ha="center", va="center")
plt.tight_layout(rect=(0.01, 0.015, .98, .98))
plt.savefig("figures/0_incidence.pdf")
plt.close()








fig, axs = plt.subplots(len(name), 1, sharex=True, figsize=(210/scalefactor, 297/scalefactor))
plt.suptitle("Seasonality", fontsize=14, y=.999)
for i in range(len(name)) :
	axs[i].plot(r[i], lw=2)
	axs[i].fill_between(range(24), rup[i], rdn[i], color=colours[0], alpha=0.3)
	axs[i].set_title(name[i], fontsize=12, loc="left")
	axs[i].set_xlim([0, 23])
	axs[i].set_xticks(range(0, 23, 4))
	axs[i].tick_params(labelsize=10)
plt.figtext(.5, 0.01, "Period", ha="center", va="center")
plt.figtext(.01, 0.5, r"Seasonality, $r$", ha="center", va="center", rotation="vertical")
plt.tight_layout(rect=(0.015, 0.015, .99, .98)) # (left, bottom, right, top) 
plt.savefig("figures/1_seasonality.pdf")
plt.close()







fig, axs = plt.subplots(len(name), 1, figsize=(210/scalefactor, 297/scalefactor))
plt.suptitle("Predicted Cases", fontsize=14, y=.999)
for i in range(len(name)) :
	axs[i].plot(t[i], ts[i], lw=2, c=colours[0])
	axs[i].plot(t[i], pred[i], lw=2, c=colours[2], alpha=0.6)
	axs[i].set_title(r"%s. Zero-corrected $R^2$ = %.02f" % (name[i], pearson[i]), fontsize=12, loc="left")
	if i == 0 :
		axs[i].legend([r"Observed Incidence", "Predicted Incidence"], loc=2)

plt.figtext(.01, 0.5, r"Incidence, $\rho C$", ha="center", va="center", rotation="vertical")
plt.figtext(.5, 0.01, "Time", ha="center", va="center")
plt.tight_layout(rect=(0.01, 0.015, .99, .98))
plt.savefig("figures/2_predictions.pdf")
plt.close()








fig, axs = plt.subplots(len(name), 1, figsize=(210/scalefactor, 297/scalefactor))
fig.suptitle("Predicted Epidemic Sizes", fontsize=14, y=.999)
for i in range(len(name)) :
	axs[i].errorbar(sizeerrx[i], sizeerry[i], yerr = sizeerre[i], fmt="o", ms=8, c=colours[0])
	axs[i].plot(sizex[i], sizey[i], lw=2, c=colours[2])
	axs[i].set_title(r"%s. $R^2$ = %.2f" % (name[i], r2[i]), fontsize=12, loc="left")

plt.figtext(.01, 0.5, "Simulated Epidemic Size", ha="center", va="center", rotation="vertical")
plt.figtext(.5, 0.01, "Actual Epidemic Size", ha="center", va="center")
plt.tight_layout(rect=(0.01, 0.015, .99, .98))
plt.savefig("figures/3_sizes.pdf")
plt.close()
