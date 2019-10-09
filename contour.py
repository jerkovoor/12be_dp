import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate

cm = plt.cm.get_cmap('RdYlBu_r')

data = np.loadtxt('/home/jerome/12Be_exp/Analysis/BeamOffset/GS_QValues_BeamOffset5_5_shift0_091_TargetDistance81_2.txt')

x = data[:, 1]
y = data[:, 2]
z = data[:, 3]

# Set up a regular grid of interpolation points
xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)

# Interpolate
rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
zi = rbf(xi, yi)

plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
           extent=[x.min(), x.max(), y.min(), y.max()], cmap=cm)
plt.colorbar()
plt.scatter(x, y, s=15, facecolors='none', edgecolors='k')
#fig, ax = plt.subplots()
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')

plt.xlabel('X [mm]')
plt.ylabel('Y [mm]')
plt.show()
