# Generate figures for paper

# Plots to generate :
#	Time-series / Predictions - alpha
#	Reporting rates
#	Periodicity
#	Predictions
#	Size

import os
import sys


sims = int(sys.argv[1])


locations = ["iceland", "bornholm", "faroe"]




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




for loc in locations :

	os.system("python batch.py %s %d" % (loc, sims))

	
	t.append(export_t)
	"""
	ts.append(export_ts)
	pred.append(export_pred)
	rho.append(export_rho)
	r.append(export_r)
	rup.append(export_rup)
	rdn.append(export_rdn)
	sizex.append(export_sizex)
	sizey.append(export_sizey)
	sizeerrx.append(export_sizeerrx)
	sizeerry.append(export_sizeerry)
	sizeerre.append(export_sizeerre)
	"""
