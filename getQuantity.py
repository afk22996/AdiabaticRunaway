import numpy as np
from read_penguin import *
import Geometry as geo
from Numerics import *

def Tau(data, planetPosition, r, N):
    azi = np.linspace(0, 2*np.pi, N)
    rGrid = cell_center(data[1])
    aziGrid = cell_center(data[2])
    
    den = np.zeros(N)
    for i in range(N):
        den[i] = findDen2D(r, azi[i], data, planetPosition)
    f = -0.5*den*np.sin(azi)
    return Simpson(f,azi[1]-azi[0])

def getTau(data, planetPosition, rLimit = (0.0, 0.3), N = 100):
    radii = np.linspace(rLimit[0], rLimit[1], N)
    tau = np.zeros(N)
    for i in range(N):
        tau[i] = Tau(data, planetPosition, radii[i], i*10+11)
    return (radii, tau)


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

def getDen2D(x, y, data, planetPosition):
    #Transform x/y to cartesian
    polar = (x,y)
    cartesian = geo.sphericalToCartesian(polar, dim = 2)
    
    #Transform to star-centric and interpolate
    starCentricCart = np.array(cartesian) + np.array(geo.sphericalToCartesian(planetPosition, dim = 2))
    starCentric = geo.cartesianToSpherical(starCentricCart, dim = 2)
    den = density(starCentric[0], starCentric[1], data)
    return den

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

def getMDot(data, rPlanet, phiPlanet, radius, N):
    angle = np.linspace(0, 2*np.pi, N)
    den = np.zeros(N)
    vrPlanet = np.zeros(N)
    
    for i in range(N):
        den[i] = findDen2D(radius,angle[i], data, (rPlanet,phiPlanet))
        vrStar,vthetaStar = getVel2D(radius, angle[i], data, (rPlanet,phiPlanet))
        cartStar = geo.sphericalToCartesian((radius, angle[i]), dim = 2) + np.array(geo.sphericalToCartesian((rplanet, phiPlanet), dim = 2))
        polarStar = np.array(geo.cartesianToSpherical(cartStar, dim = 2))
        vrPlanet[i] = getvr2D(vrStar, vthetaStar,polarStar[0], polarStar[1], rPlanet, phiPlanet)
    f = radius*den*vrPlanet
    return Simpson(f,angle[1]-angle[0])

def findFg(i, j, data, planetPosition, planetMass):
    rf = data[1][i+1]
    ri = data[1][i]
    r = 0.5*(rf + ri)
    
    phif = data[2][j+1]
    phii = data[2][j]
    phi = 0.5*(phif+phii)
    
    den = data[3][j,i]
    
    area = (phif-phii)*(rf**2 - ri**2)/(2)
    mass = den*area
    
    cart = np.array(geo.sphericalToCartesian((r, phi), dim = 2))
    cartPlanet = cart - np.array(geo.sphericalToCartesian(planetPosition, dim = 2))
    polarPlanet = np.array(geo.cartesianToSpherical(cartPlanet, dim = 2))
    
    Fg = -(mass*planetMass)/(polarPlanet[0]**2)
    return np.array([Fg*np.cos(phi), Fg*np.sin(phi)])

def findTau(i, j, data, planetPosition, planetMass):
    FgPhi = findFg(i, j, data, planetPosition, planetMass)[1]
    tau = -FgPhi*planetPosition[0]
    return tau

def getTau(data, planetPosition, planetMass):
    tau = 0
    for i in range(len(data[1]) - 1):
        for j in range(len(data[2]) - 1):
            tau = tau + findTau(i, j, data, planetPosition, planetMass)
    return tau

def getFg(data, planetPosition, planetMass):
    Fg = np.array([0,0])
    for i in range(len(data[1]) - 1):
        for j in range(len(data[2]) - 1):
            Fg = Fg + findFg(i, j, data, planetPosition, planetMass)
    return Fg