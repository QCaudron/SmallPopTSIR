import os
import pandas as pd
import numpy as np
import seaborn

c = seaborn.color_palette("deep", 8)

os.chdir("../")

sumR = []

for i in range(1, 25) :

	os.system("python batch.py iceland 500 %d" % i)
	os.system("python batch.py faroe 500 %d" % i)
	os.system("python batch.py bornholm 500 %d" % i)

	r = []

	for file in os.listdir("paper/figures/") :
		if file.endswith(".json") :
			data = pd.read_json("paper/figures/%s" % file)
			r.append(data["pearson"].values[0])


	print "-------- %d" % i

	sumR.append(np.sum(r))


print "BEST SENSITIVITY : %d" % (np.argmax(sumR) + 1)




rvik = [-21.9333, 64.1333]
akur = [-18.1000, 65.6833]
vest = [-20.2833, 63.4167]
haf = [-21.9542, 64.0678]
m = Basemap(projection="merc", llcrnrlat=63.1, urcrnrlat=66.7, llcrnrlon=-25.1, urcrnrlon=-12.9, resolution="f");
m.drawcoastlines(linewidth=2)
m.drawrivers()
m.shadedrelief()

x, y = m(rvik[0], rvik[1])
m.scatter(x, y, 80, marker="o", color=c[2])

x, y = m(akur[0], akur[1])
m.scatter(x, y, 80, marker="o", color=c[2])

x, y = m(vest[0], vest[1])
m.scatter(x, y, 80, marker="o", color=c[2])

x, y = m(haf[0], haf[1])
m.scatter(x, y, 80, marker="o", color=c[2])

plt.savefig("/Users/qcaudron/SmallPopTSIR/talk/iceland.pdf")



m = Basemap(projection="merc", llcrnrlat=54.94, urcrnrlat=55.33, llcrnrlon=14.65, urcrnrlon=15.2, resolution="f");
m.drawcoastlines(linewidth=2)
m.drawrivers()
m.shadedrelief()

plt.savefig("/Users/qcaudron/SmallPopTSIR/talk/bornholm.pdf")



m = Basemap(projection="merc", llcrnrlat=61.336, urcrnrlat=62.423545, llcrnrlon=-7.8, urcrnrlon=-6.099557, resolution="f");
m.drawcoastlines(linewidth=2)
m.drawrivers()
m.shadedrelief()

plt.savefig("/Users/qcaudron/SmallPopTSIR/talk/faroe.pdf")




m = Basemap(projection="merc", llcrnrlat=66.912917, urcrnrlat=53.304748, llcrnrlon=-25.743112, urcrnrlon=17.762748, resolution="f");
m.drawcoastlines(linewidth=2)
m.drawrivers()
m.bluemarble()
plt.savefig("/Users/qcaudron/SmallPopTSIR/talk/big.pdf")
