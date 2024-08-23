import numpy as np
import Geometry as geo
from Numerics import RK5Step, binSearch

def findH(r, fun, coordX, coordY, coordZ, direction):
    x = r[0]
    y = r[1]
    z = r[2]
    
    vx, vy, vz = fun(x, r)
    xPoints = binSearch(coordX, 0, len(coordX), x)
    yPoints = binSearch(coordY, 0, len(coordY), y)
    zPoints = binSearch(coordZ, 0, len(coordZ), z)
    if(xPoints[0] == -np.infty):
        return 0
    elif(xPoints[1] == np.infty):
        return 0
    if(yPoints[0] == -np.infty):
        yPoints = (len(coordY)-1, 0)
    elif(yPoints[1] == np.infty):
        yPoints = (0, len(coordY)-1)
    if(zPoints[1] == np.infty):
        if z > np.pi/2:
            z = np.pi - z
            zPoints = binSearch(coordZ, 0, len(coordZ), z)
            if(zPoints[0] == -np.infty):
                return 0
            if(zPoints[1] == np.infty):
                return 0
        else:
            zPoints = (-1, -1)
    if(zPoints[0] == -np.infty):
        return 0
    xp = xPoints[direction]
    yp = yPoints[direction]
    zp = zPoints[direction]
    if(vx != 0):
        hx = abs(coordX[xp] - x)/vx
    else:
        hx = np.infty
    if(vy != 0):
        hy = abs(coordY[yp] - y)/vy
    else:
        hy = np.infty
    if(vz != 0):
        hz = abs(coordZ[zp] - z)/vz
    else:
        hz = np.infty
    if(hx == 0 and hy == 0 and hz == 0):
        return 0
    if(hx == 0):
        hx = np.infty
    if(hy == 0):
        hy = np.infty
    if(hz == 0):
        hz = np.infty
    h = min(abs(hx), abs(hy), abs(hz))
    return h
    

'''
Algorithm to return the flow line of a particle in a 2d velocity vector field

Inputs:
xi - x point which the flow line passes through
yi - y point which the flow line passes through
Xs - Collection of grid points in the X direction
Ys - Collection of grid points in the Y direction
maxerror - Desired error for the flow lines
fun - function which returns the x and y velocities at a given point

Returns:
xs - x points which the flow line passes through
ys - y points which the flow line passes through
'''
def flowLine2D(xi, yi, Xs, Ys, maxerror, fun):
    ys = []
    xs = []
    y0 = np.array([xi, yi])
    xf = np.max(Xs)
    yf = np.max(Ys)

    h = abs(xf-xi)
    ys.append(yi)
    xs.append(xi)
    n = 0
    aveCellX = (xf - xi)/len(Xs)
    aveCellY = (yf - yi)/len(Ys)
    while(abs(y0[0]) < xf and abs(y0[1]) < yf):
        correction = RK5Step(y0, y0[0], xf, h, maxerror, fun)
        if(abs(correction[0][0]) <= 0.01*aveCellX and abs(correction[0][1]) <= 0.01*aveCellY):
            break
        elif(n > 1000):
            break
        y0 = y0 + correction[0]
        ys.append(y0[1])
        xs.append(y0[0])
        h = correction[3]
        n = n+1
    y0 = np.array([xi, yi])
    h = abs(xf-xi)
    n = 0
    while(abs(y0[0]) < xf and abs(y0[1]) < yf):
        correction = RK5Step(y0, y0[0], xf, h, maxerror, fun)
        if(abs(correction[0][0]) <= 0.01*aveCellX and abs(correction[0][1]) <= 0.01*aveCellY):
            break
        elif(n > 1000):
            break
        y0 = y0 - correction[0]
        ys.insert(0, y0[1])
        xs.insert(0, y0[0])
        h = correction[3]
        n = n+1
    return (xs, ys)


'''
Algorithm to return the flow line of a particle in a 3d velocity vector field

Inputs:
xi - x point which the flow line passes through
yi - y point which the flow line passes through
zi - z point which the flow line passes through
Xs - Collection of grid points in the X direction
Ys - Collection of grid points in the Y direction
Zs - Collection of grid points in the Z direction
maxerror - Desired error for the flow lines
fun - function which returns the x and y velocities at a given point

Returns:
xs - x points which the flow line passes through
ys - y points which the flow line passes through
zs - z points which the flow line passes through
'''
def flowLine3D(xi, yi, zi, Xs, Ys, Zs, maxerror, fun, maxstep):
    zs = []
    xs = []
    ys = []
    y0 = np.array([xi, yi, zi])
    xf = np.max(np.abs(Xs))
    yf = np.max(np.abs(Ys))
    zf = np.max(np.abs(Zs))
    zs.append(zi)
    xs.append(xi)
    ys.append(yi)
    n = 0
    phiInitial = geo.cartesianToSpherical(y0, dim = 3)[1]
    lastPhi = phiInitial
    planetCoords = (1, np.pi, np.pi/2)
    h = findH(y0, fun, Xs, Ys, Zs, 1)
    while(abs(y0[0]) <= abs(xf) and abs(y0[1]) <= abs(yf) and abs(y0[2]) <= abs(zf)):
        correction = RK5Step(y0, y0[0], xf, h, maxerror, fun)
        y0 = y0 + correction[0]
        phi = geo.cartesianToSpherical(y0, dim = 3)[1]
        zs.append(y0[2])
        xs.append(y0[0])
        ys.append(y0[1])
        h = findH(y0, fun, Xs, Ys, Zs, 1)
        if(h == 0 or h > 5):
            #print(h, y0)
            break
        phi = geo.cartesianToSpherical(y0, dim = 3)[1]
        if(n > 0):
            if((phi <= phiInitial and lastPhi >= phiInitial) or (phi >= phiInitial and lastPhi <= phiInitial)):
                #print(phi, phiInitial)
                break
        lastPhi = phi
        n = n+1

    '''y0 = np.array([xi, yi, zi])
    h = findH(y0, fun, Xs, Ys, Zs, 0)
    n = 0
    lastPhi = phiInitial
    h = findH(y0, fun, Xs, Ys, Zs, 0)
    while(abs(y0[0]) <= abs(xf) or abs(y0[1]) <= abs(yf) or abs(y0[2]) <= abs(zf)):
        correction = RK5Step(y0, y0[0], xf, h, maxerror, fun)
        if(h == 0):
            print(h, y0)
            break
        y0 = y0 - correction[0]
        zs.insert(0, y0[2])
        xs.insert(0, y0[0])
        ys.insert(0, y0[1])
        phi = geo.cartesianToSpherical(y0, dim = 3)[1]
        if(n > 0):
            if(n > 0):
                if((phi <= phiInitial and lastPhi >= phiInitial) or (phi >= phiInitial and lastPhi <= phiInitial)):
                    print(y0)
                    break
        lastPhi = phi
        h = findH(y0, fun, Xs, Ys, Zs, 0)
        n = n+1'''
    return (xs, ys, zs)

def isHorseshoe(flow, planetCoords):
    xc,yc,zc = geo.sphericalToCartesian(planetCoords, dim = 3)
    if(len(flow[0]) <= 1):
        return False
    radii = [np.sqrt((flow[0][i] + xc)**2 + (flow[1][i] + yc)**2) for i in range(int(len(flow[0])))]
    if(np.min(radii) <= planetCoords[0] and np.max(radii) >= planetCoords[0]):
        return True
    else:
        return False