import numpy as np
import Geometry as geo
from Interpolate import *
from Search import binSearch
from read_penguin import *
from mass_flux_2D import *

'''Algorithm to return the mean temperature for a circle or radius r around the planet. Intended to be a helper function to findTemp, which returns the temperatures along a given range but can be used individually

inputs:
dat - simulation data
planet_r - radial position of the planet
planet_phi - angular position of the planet
rad - radial distance from the planet
N - Number of points to be considered for the mean'''

def getTemp(dat,planet_r,planet_phi,rad,N):

    azi = np.linspace(0.0,2.0*np.pi,num=N)
    r_star, phi_star = planet_to_star(rad,azi,planet_r,planet_phi)

    r_grid = cell_center(dat[1])
    phi_grid = cell_center(dat[2])

    temp = np.zeros(N)
    for i in range(N):
        density = interpolate_2D_periodicY((r_star[i],phi_star[i]),r_grid,phi_grid,dat[3],dat[2][0],dat[2][-1])
        pressure = interpolate_2D_periodicY((r_star[i],phi_star[i]),r_grid,phi_grid,dat[4],dat[2][0],dat[2][-1])
        temp[i] = pressure/density
    return np.mean(temp)

'''Algorithm to return the mean temperatures of increasing radii around the planet

inputs:
data - simulation data
rPlanet - radial position of the planet
phiPlanet - angular position of the planet
rlim - interval over which the temperature should be evaluated
N - number of desired points in the interval

returns:
r - radii over which the planet was evaluated
T - average temperature at these corresponding radii
'''
def findTemp(data, rPlanet, phiPlanet, rlim = (0, 0.3), N = 100):
    r = np.linspace(rlim[0],rlim[1],N)
    T = np.zeros(N)
    for i in range(N):
        T[i] = getTemp(data, rPlanet, phiPlanet, r[i], i*10+11)
    return (r, T)

'''Algorithm to interpolate 2D density data with polar boundary conditions (old code but still used in getDen2D so I haven't removed it, do not use)'''
def density(x, y, data):
    #Parsing relevant parts of data input
    xVals = cell_center(data[1])
    yVals = cell_center(data[2])
    dens = data[3]
    
    #Finding indices in the data between which the target points lie
    xPoints = binSearch(xVals, 0, len(xVals), x)
    yPoints = binSearch(yVals, 0, len(yVals), y)
    
    #Applying Polar Boundary Conditions (y in this case is angle, so it should be periodic but x is radius so anything outside of data set is made to be 0)
    if(xPoints[0] == -np.infty):
        return 0
    elif(xPoints[1] == np.infty):
        return 0
    if(yPoints[0] == -np.infty):
        yPoints = (0, len(yVals)-2)
    elif(yPoints[1] == np.infty):
        yPoints = (len(yVals)-2, 0)
    
    #Setting up grid points for interpolation (these are the indices of the square lattice points)
    lowx = xPoints[0]
    highx = xPoints[1]
    lowy = yPoints[0]
    highy = yPoints[1]
    
    #Creating values at each square lattice point and interpolating
    squareVals = [dens[lowy,lowx], dens[highy,lowx], dens[lowy,highx], dens[highy,highx]]
    targetCoords = (x,y)
    minCoords = (xVals[lowx], yVals[lowy])
    maxCoords = (xVals[highx], yVals[highy])
    den = biInterpolate(targetCoords, squareVals, minCoords, maxCoords)
    return den

'''Algorithm to find density at a given point (assumed planet-centric) in a star-centric grid

inputs:
x - radius in planet-centric grid
y - angle in planet-centric grid
data - simulation data
planetPosition - indexable object containing the radial and angular positions of the planet

returns:
den - interpolated density at the desired point'''
def getDen2D(x, y, data, planetPosition):
    #Transform x/y to cartesian
    polar = (x,y)
    cartesian = geo.sphericalToCartesian(polar, dim = 2)
    
    #Transform to star-centric and interpolate
    starCentricCart = np.array(cartesian) + np.array(geo.sphericalToCartesian(planetPosition, dim = 2))
    starCentric = geo.cartesianToSpherical(starCentricCart, dim = 2)
    den = density(starCentric[0], starCentric[1], data)
    return den

'''Algorithm to find and transorm the velocity at a given point (assumed planet-centric) in a 2D polar grid

inputs:
x - radius in planet-centric grid
y - angle in planet-centric grid
data - simulation data
planetPosition - indexable object containing the radial and angular positions of the planet
planetVel - indexable object containing the radial and angular velocities of the planet
corot - boolean indicating if you want the velocities in the corotating frame (defaults to true)
cart - boolean indicating if you want to velocities to be outputted as cartesian values instead of polar (defaults to false)

returns:
v - array containing the interpolated 2D velocity data (either polar or cartesian)'''
def getVel2D(x, y, data, planetPosition, planetVel, corot = True, cart = False):
    rGrid = cell_center(data[1])
    phiGrid = cell_center(data[2])
    
    planetPolar = (x,y)
    planetCart = geo.sphericalToCartesian(planetPolar, dim = 2)
    starCart = np.array(planetCart) + np.array(geo.sphericalToCartesian(planetPosition, dim = 2))
    starPolar = geo.cartesianToSpherical(starCart, dim = 2)
    
    vrStar = interpolate_2D_periodicY(starPolar, rGrid, phiGrid, data[5], data[2][0], data[2][-1])
    vphiStar = interpolate_2D_periodicY(starPolar, rGrid, phiGrid, data[6], data[2][0], data[2][-1])
    if(vrStar == None or vphiStar == None):
        return [0,0]
    if(corot):
        vphiStar = vphiStar - starPolar[0]*planetVel[1]
    v = np.array([vrStar, vphiStar])
    if(cart):
        return geo.sphericalToCartesianVelocity(starPolar, v, dim = 2)
    else:
        return v

'''Algorithm to find and transorm the velocity at a given point (assumed planet-centric) in a 3D spherical grid

inputs:
x - radius in planet-centric grid
y - polar angle in planet-centric grid
z - azimuthal angle in planet-centric grid
data - simulation data
planetPosition - indexable object containing the radial and angular positions of the planet
planetVel - indexable object containing the radial and angular velocities of the planet
corot - boolean indicating if you want the velocities in the corotating frame (defaults to true)
cart - boolean indicating if you want to velocities to be outputted as cartesian values instead of spherical (defaults to false)

returns:
v - array containing the interpolated 3D velocity data (either spherical or cartesian)'''
def getVel3D(coords, planetCoords, planetVel, data, corot = True, cart = False):
    xVals = cell_center(data[1])
    yVals = cell_center(data[2])
    zVals = cell_center(data[3])
    xVel = data[6]
    yVel = data[7]
    zVel = data[8]
    planetSphere = coords
    planetCart = geo.sphericalToCartesian(planetSphere, dim = 3)
    starCart = np.array(planetCart) + np.array(geo.sphericalToCartesian(planetCoords, dim = 3))
    starSphere = geo.cartesianToSpherical(starCart, dim = 3)
    vx = interpolate3DSpherical(xVals,yVals,zVals,xVel,starSphere)
    vy = interpolate3DSpherical(xVals,yVals,zVals,yVel,starSphere)
    if(corot):
        vy = vy - starSphere[0]*np.sin(starSphere[2])
    vz = np.sign(np.pi/2 - starSphere[2])*interpolate3DSpherical(xVals,yVals,zVals,zVel,starSphere)
    v = np.array([vx, vy, vz])
    if(cart):
        return np.array(geo.sphericalToCartesianVelocity(starSphere, v, dim = 3))
    else:
        return v

'''Algorithm to transform velocity in a star-centric frame to radial velocity in a planet-centric frame

inputs:
vr_star - star-centric radial velocity
vp_star - star-centric angular velocity
r_star - star-centric radial position
phi_star - star-centric angular position
r_planet - radial position of the planet in star-centric frame
phi_planet - angular position of the planet in star-centric frame

returns:
vr_planet - radial velocity in planet's reference frame'''

def getvr2D(vr_star,vp_star,r_star,phi_star,r_planet,phi_planet):
    polarvel = [vr_star, vp_star]
    cartvel = geo.sphericalToCartesianVelocity((r_star, phi_star), polarvel, dim = 2)

    vr_planet = cartvel[0]*np.cos(phi_planet) + cartvel[1]*np.sin(phi_planet)
    return vr_planet

'''Algorithm to transform velocity in a star-centric frame to radial velocity in a planet-centric frame (Note: Untested as of writing)

inputs:
vr_star - star-centric radial velocity
vp_star - star-centric azimuthal velocity
vt_star - star-centric polar velocity
r_star - star-centric radial position
phi_star - star-centric azimuthal position
theta_star - star-centric polar position
r_planet - radial position of the planet in star-centric frame
phi_planet - azimuthal position of the planet in star-centric frame
theta_planet - polar position of the planet in star-centric frame

returns:
vr_planet - radial velocity in planet's reference frame'''
def getvr3D(vr_star, vp_star, vt_star, r_star, phi_star, theta_star, r_planet, phi_planet, theta_planet):
    spherical = [r_star, phi_star, theta_star]
    sphereVel = [vr_star, vphi_star, vt_star]
    cartVel = geo.sphericalToCartesianVelocity(spherical, sphereVel, dim = 3)
    
    vr_planet = cartVel[0]*np.cos(phi_planet)*np.sin(theta_planet) + cartVel[1]*np.sin(phi_planet)*np.sin(theta_planet) + cartVel[2]*np.cos(theta_planet)
    return vr_planet