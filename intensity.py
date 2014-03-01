import numpy as np
import matplotlib.pyplot as pl
import scipy.interpolate

x, y, z = np.loadtxt("data_temperaturemap.txt").transpose()

# Set up a regular grid of interpolation points
xi, yi = np.linspace(x.min(), x.max(), 300), np.linspace(y.min(), y.max(), 300)
xi, yi = np.meshgrid(xi, yi)

# Interpolate; there's also method='cubic' for 2-D data such as here
zi = scipy.interpolate.griddata((x, y), z, (xi, yi), method='linear')

pl.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
           extent=[x.min(), x.max(), y.min(), y.max()])
pl.scatter(x, y, c=z)
pl.colorbar()

x, y, vx, vy = np.loadtxt("data_flux.txt").transpose()
pl.quiver(x,y,vx,vy, width=0.002)
pl.show()