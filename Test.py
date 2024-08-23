import Geometry as geo
from getQuantity import getVel3D
from Flow import flowLine3D, isHorseshoe, findH
import numpy as np
from read_penguin import load_3D_data, cell_center
from matplotlib import pyplot as plt
from time import time
from Numerics import *

fpath = "/home/alex/Documents/Coding/Research/AdiabaticRunaway/Data/"
xres = 288
yres = 480
zres = 144
label = "h50_1p10E_e0_PPM4_ave"
frame = 10
data3d = load_3D_data(fpath, xres, yres, zres, label, frame)

planetCoords = (1,np.pi,np.pi/2)
xp, yp, zp = geo.sphericalToCartesian(planetCoords, dim = 3)
planetVel = (0,1,0)

planetMass = 0.00003
h = 0.05
gamma = 1.4

error = 1e-8
maxstep = 0.01

def isoVel3D(x, y):
    coords = geo.cartesianToSpherical(y, dim = 3)
    v = getVel3D(coords, planetCoords, planetVel, data3d, cart = True)
    return v

rGrid = cell_center(data3d[1])
phiGrid = cell_center(data3d[2])
thetaGrid = cell_center(data3d[3])
coordX = np.linspace(-rGrid[-1] - xp, rGrid[-1] - xp, 1000)
coordY = np.linspace(-rGrid[-1] - yp, rGrid[-1] - yp, 1000)
coordZ = np.linspace(-np.cos(thetaGrid[0]) - zp, np.cos(thetaGrid[0]) - zp, 1000)

def starToPlanet(sCoords, planetCoords):
    sCoordsCart = geo.sphericalToCartesian(sCoords, dim = 3)
    pCoordsCart = np.array(sCoordsCart) - np.array(geo.sphericalToCartesian(planetCoords, dim = 3))
    return pCoordsCart
def velMag(r):
    phi = 2*np.pi/3
    theta = np.pi/2
    x,y,z = starToPlanet((r,phi,theta), planetCoords)
    vx,vy,vz = isoVel3D(0, (x,y,z))
    return np.sqrt(vx**2 + vy**2 + vz**2)


corot = 0.9969806020039222
guess = 0.01
xi, yi, zi = starToPlanet((corot - guess, 2*np.pi/3, thetaGrid[-2]), planetCoords)
flow = flowLine3D(xi, yi, zi, coordX, coordY, coordZ, error, isoVel3D, maxstep)

ax = plt.figure().add_subplot(projection = '3d')
ax.plot(flow[0], flow[1], flow[2])
ax.set_zlim(0, 0.1)
plt.show()