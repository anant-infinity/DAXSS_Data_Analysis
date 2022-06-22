from xspec import *
import matplotlib.pyplot as plt

# Clearing Old Data + Models
AllData.clear()
AllModels.clear()

# Setting Feldman Abundances
#Xset.abund = 'feld'
Xset.show()
# Loading the DAXSS spectrum
spec = Spectrum("/home/anant/PycharmProjects/DAXSS_Data_Analysis/minxss_fm3_PHA_2022-03-31T19-07-48Z.pha")
spec.ignore('**-0.7 3.5-**')

# define the model
m1 = Model("vvapec")

# Free some parameters that are frozen (Mg, Si, and S)

m1.vvapec.Mg.frozen=False
m1.vvapec.Si.frozen=False
m1.vvapec.Fe.frozen=False
m1.vvapec.S.frozen=False
m1.vvapec.Ne.frozen=False
#m1.vvapec.Ca.frozen=False

# do the fit
Fit.nIterations = 100
Fit.perform()
# plot data, model and del-chi
# Plot.device = '/xw'
Plot.xAxis = 'keV'
Plot('ld','delc')
ene = Plot.x(plotGroup=1, plotWindow=1)
eneErr = Plot.xErr(plotGroup=1, plotWindow=1)
spec = Plot.y(plotGroup=1, plotWindow=1)
specErr = Plot.yErr(plotGroup=1, plotWindow=1)

fitmodel = Plot.model(plotGroup=1, plotWindow=1)

delchi = Plot.y(plotGroup=1, plotWindow=2)
delchiErr = Plot.yErr(plotGroup=1, plotWindow=2)

fig0 = plt.figure(num=None, figsize=(6, 4), facecolor='w', edgecolor='k')

ax0 = fig0.add_axes([0.15, 0.4, 0.8, 0.55])
ax0.xaxis.set_visible(False)
plt.errorbar(ene, spec, xerr=eneErr, yerr=specErr, fmt='.', ms=0.5, capsize=1.0, lw=0.8,label = 'Data')
plt.step(ene, fitmodel, where='mid',label = 'Model')
plt.yscale("log")
plt.xlim([0.7, 3.5])
#plt.ylim([1, 1e6])
plt.legend()
plt.ylabel('Rate (counts s$^{-1}$ keV$^{-1}$)')

ax1 = fig0.add_axes([0.15, 0.15, 0.8, 0.25])
plt.axhline(0, linestyle='dashed', color='black')
plt.errorbar(ene, delchi, xerr=eneErr, yerr=delchiErr, fmt='.', ms=0.1, capsize=1.0, lw=0.8)
plt.xlim([0.7, 3.5])
plt.ylabel('$\Delta \chi$')
plt.xlabel('Energy (keV)')
plt.show()
plt.close()