from read_penguin import *
import matplotlib.pyplot as plt
import numpy as np
from Transform import *
import multiprocessing as mp
from mass_flux_2D import *
from Flow import *
import Geometry as geo
from Search import *
from time import time
from getQuantity import *
from Interpolate import *

def isoVel3D(x, y):
    planetCoords = (1,0,np.pi/2)
    planetVel = (0,1,0)
    coords = geo.cartesianToSpherical(y, dim = 3)
    return np.array([0,0, getVel3D(coords, planetCoords, planetVel, data3d, cart = True)[2]])

def isHorseshoe(flow, planetCoords):
    xc,yc,zc = geo.sphericalToCartesian(planetCoords, dim = 3)
    if(len(flow[0]) <= 1):
        return False
    radii = [np.sqrt((flow[0][i] - xc)**2 + (flow[1][i] - yc)**2) for i in range(int(0.9*len(flow[0])))]
    if(np.min(radii) < planetCoords[0] and np.max(radii) > planetCoords[0]):
        return True
    else:
        return False

binpath = "C:\\Users\\sidio\\OneDrive\\Desktop\\Coding\\Analysis\\"
xres = 480
yres = 864
zres = 264
label = "h50_1p10E_e0_PPM4_ave"


data = load_3D_data(binpath, xres, yres, zres, label, 0)

planetCoords = (1, 0, np.pi/2)
planetVel = (0,1,0)

rGrid = cell_center(data[1])
phiGrid = cell_center(data[2])
thetaGrid = cell_center(data[3])

phiPoints = np.linspace(0, np.pi, 1000)
dens = [interpolate3DSpherical(rGrid, phiGrid, thetaGrid, data[4], (1, 0, phiPoints[i])) for i in range(len(phiPoints))]

plt.plot(phiPoints, dens)
plt.xticks([0, np.pi/2, np.pi])
plt.vlines(np.pi/2, 0, 1, linestyles='dashed')
plt.show()