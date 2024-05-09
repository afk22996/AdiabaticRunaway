import numpy as np
import Geometry as geo
from Interpolate import *
from Search import binSearch
from read_penguin import *

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

def findDen2D(x, y, data, planetPosition):
    #Transform x/y to cartesian
    polar = (x,y)
    cartesian = geo.sphericalToCartesian(polar, dim = 2)
    
    #Transform to star-centric and interpolate
    starCentricCart = np.array(cartesian) + np.array(geo.sphericalToCartesian(planetPosition, dim = 2))
    starCentric = geo.cartesianToSpherical(starCentricCart, dim = 2)
    den = density(starCentric[0], starCentric[1], data)
    return den

def findVel2D(x, y, data, planetPosition, planetVel, cart = True):
    if(cart):
        cartPlanet = np.array([x,y])
        cartStar = cartPlanet + np.array(geo.sphericalToCartesian(planetPosition, dim = 2))
        polarStar = geo.cartesianToSpherical(cartStar, dim = 2)
        if(polarStar[1] > 2*np.pi):
            polarStar[1] = polarStar[1]%(2*np.pi)
        xPoints = binSearch(data[1], 0, len(data[1]), polarStar[0])
        yPoints = binSearch(data[2], 0, len(data[2]), polarStar[1])
        if(xPoints[0] == -np.infty or xPoints[1] == np.infty):
            return(0,0)
        if(yPoints[0] == -np.infty):
            yPoints = (0, len(data[2])-2)
        if(yPoints[1] == np.infty):
            yPoints = (len(data[2])-2, 0)
        lowx = xPoints[0]
        lowy = yPoints[0]
        highx = xPoints[1]
        highy = yPoints[1]

        xSquareVals = [data[5][lowy,lowx], data[5][highy,lowx], data[5][lowy, highx], data[5][highy, highx]]
        ySquareVals = [data[6][lowy,lowx], data[6][highy,lowx], data[6][lowy, highx], data[6][highy, highx]]

        targetCoords = (polarStar[0],polarStar[1])
        minCoords = (data[1][lowx], data[2][lowy])
        maxCoords = (data[1][highx], data[2][highy])
        vx = biInterpolate(targetCoords, xSquareVals, minCoords, maxCoords)
        vy = biInterpolate(targetCoords, ySquareVals, minCoords, maxCoords)
        vStar = np.array([vx, vy])
        vStar[1] = vStar[1] - planetVel[1]*polarStar[0]
        return geo.sphericalToCartesianVelocity(polarStar, vStar, dim = 2)
    
    else:
        cartPlanet = geo.sphericalToCartesian((x,y), dim = 2)
        xPoints = binSearch(data[0], 0, len(data[0]), cartPlanet[0])
        yPoints = binSearch(data[1], 0, len(data[1]), cartPlanet[1])
        if(xPoints[0] == -np.infty or xPoints[1] == np.infty):
            return(0,0)
        if(yPoints[0] == -np.infty or yPoints[1] == np.infty):
            return(0,0)
        lowx = xPoints[0]
        lowy = yPoints[0]
        highx = xPoints[1]
        highy = yPoints[1]

        xSquareVals = [data[2][lowy,lowx], data[2][highy,lowx], data[2][lowy, highx], data[2][highy, highx]]
        ySquareVals = [data[3][lowy,lowx], data[3][highy,lowx], data[3][lowy, highx], data[3][highy, highx]]

        targetCoords = (cartPlanet[0],cartPlanet[1])
        minCoords = (data[0][lowx], data[1][lowy])
        maxCoords = (data[0][highx], data[1][highy])
        vx = biInterpolate(targetCoords, xSquareVals, minCoords, maxCoords)
        vy = biInterpolate(targetCoords, ySquareVals, minCoords, maxCoords)
        vStar = np.array(geo.cartesianToSphericalVelocity(cartPlanet, [vx, vy], dim = 2))
        return vStar
    
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

def findTemp(data, rPlanet, phiPlanet):
    r = np.linspace(0,0.3,100)
    T = np.zeros(100)
    for i in range(100):
        T[i] = getTemp(data, rPlanet, phiPlanet, r[i], i*10+11)
    return (r, T)