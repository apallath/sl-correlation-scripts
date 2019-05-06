# Code to plot mean square displacements output from LAMMPS

import numpy as np
import matplotlib.pyplot as plt

# REQUIRED INPUT
kvals = [0.001,0.1,1,10,100,1000,10000]

msdmeanplot = []
msdstdplot = []

for k in kvals:
    f = open("MSD_all_{}.dat".format(k), "r")
    f.readline()
    msd = []
    for line in f:
        [msdx,msdy,msdz,msdw] = map(float,line.strip().split())
        msd.append(msdw)
    msdmean = np.mean(msd)
    msdstd = np.std(msd)
    msdmeanplot.append(msdmean)
    msdstdplot.append(msdstd)

plt.figure()
plt.errorbar(kvals,msdmeanplot,msdstdplot)
ax = plt.gca()

#set log scale
ax.set_xscale('log')
ax.set_yscale('log')
plt.title("LJ potential, $\sigma = 1$, $\epsilon = 0.05$, $r_{cut} = 2.5$")
plt.xlabel("Spring constant (F/$\sigma$)")
plt.ylabel("Mean square displacement ($\sigma$)")

#mark points
for i, txt in enumerate(msdmeanplot):
    ax.annotate("{:.2E}".format(txt), (kvals[i], txt))  

# OUTPUTS
plt.savefig("MSD.pdf")
plt.savefig("MSD.png", dpi=120)
