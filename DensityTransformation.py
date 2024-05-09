'''
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
This file overall transforms densities from a star-centric 2D polar grid to a star-centric 2D cartesian grid, then to a planet-centric 2D polar grid. It is included on this github as an example of how to transform an invariant value, to look at
transforming a value which needs to be transformed itself (i.e. not just change of coordinates) such as velocity, the Jupyter notebook file does that for velocities
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
'''
from read_penguin import load_2D_data
import matplotlib.pyplot as plt
import numpy as np
import geometry as geo
from Interpolate import biInterpolate
from Search import binSearch

'''
Function to interpolate the densities from a given cartesian coordinate from a polar data source

Inputs:
(x,y) - coordinates of interest
data - read_penguin output (assumes 2D, would need to be adapted for 3D use)
'''
def density(x, y, data):
    #Parsing relevant parts of data input
    xVals = data[1]
    yVals = data[2]
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
    
'''
Very similar to density, but assumes cartesiaan boundary conditions instead of polar (i.e. everything outside of data range is 0). The different inputs are due to not assuming that the data coming is from read_penguin, now you just give
the x/y indices and the data points. The purpose of this is to interpolate from the 2D cartesian grid to another 2D polar grid
'''
def cartDensity(x, y, xVals, yVals, dens):
    xPoints = binSearch(xVals, 0, len(xVals), x)
    yPoints = binSearch(yVals, 0, len(yVals), y)
    if(xPoints[0] == -np.infty or yPoints[0] == -np.infty):
        return 0
    if(xPoints[1] == np.infty or yPoints[1] == np.infty):
        return 0
    lowx = xPoints[0]
    highx = xPoints[1]
    lowy = yPoints[0]
    highy = yPoints[1]
    squareVals = [dens[lowy,lowx], dens[highy,lowx], dens[lowy,highx], dens[highy,highx]]
    targetCoords = (x,y)
    minCoords = (xVals[lowx], yVals[lowy])
    maxCoords = (xVals[highx], yVals[highy])
    den = biInterpolate(targetCoords, squareVals, minCoords, maxCoords)
    return den

if __name__ == "__main__":
    #Loading density data
    xres = 384
    yres = 768
    filepath = "/home/afkirby/penguinPlots/2DGammaChange/"
    gam10data = load_2D_data("/scratch/afkirby/2DAdiabaticParameterChange/Gamma1.0/", xres, yres, "h50_1p1J_e0_PPM4", 100)

    #Creating cartesian grid/defining planet position in both coordinate systems
    coordX = np.ndarray(1001)
    coordY = np.ndarray(1001)
    planetCoords = (1, np.pi)
    planetCoordsCart = geo.sphericalToCartesian(planetCoords, dim = 2)

    #Populating Cartesian Grid
    for i in range(1001):
        coordX[i] = -gam10data[1][-1] + 2*gam10data[1][-1]*(i)/1000

    for j in range(1001):
        coordY[j] = -gam10data[1][-1] + 2*gam10data[1][-1]*(j)/1000

    #Interpolating over polar array to find data for cartesian coordinates
    gam10denCart = np.ndarray((1001,1001))
    for i in range(1001):
        for j in range(1001):
            cartesian = (coordX[i], coordY[j])
            polar = geo.cartesianToSpherical(cartesian, 2)
            gam10denCart[j,i] = density(polar[0], polar[1], gam10data)

    #Shifting Cartesian Grid to be planet-centric
    for i in range(1001):
        coordX[i] = coordX[i] - planetCoordsCart[0]
        coordY[i] = coordY[i] - planetCoordsCart[1]

    #Creating planet-centric polar grid
    gam10denplanet = np.ndarray((xres, yres))
    planetR = np.ndarray(xres+1)
    planetTheta = np.ndarray(yres+1)
    for i in range(xres+1):
        planetR[i] = 2.7*i/(xres)
    for j in range(yres+1):
        planetTheta[j] = 2*np.pi*j/(yres)

    #Interpolating over cartesian grid to find data for the new, planet-centric, polar grid
    for i in range(xres):
        for j in range(yres):
            polar = (planetR[i], planetTheta[j])
            cartesian = geo.sphericalToCartesian(polar, 2)
            gam10denplanet[i,j] = cartDensity(cartesian[0], cartesian[1], coordX, coordY, gam10denCart)

    #Plotting

    plt.figure()
    plt.pcolor(planetTheta, planetR, gam10denplanet)
    plt.title("Final Planet-Centric Isothermal Density")
    plt.xlabel("Theta (Rad)")
    plt.ylabel("R (Planet a)")
    plt.colorbar()
    plt.savefig(filepath + "Gamma1.0FinalDensityPolarPlanet.png")