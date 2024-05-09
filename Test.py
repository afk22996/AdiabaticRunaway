from read_penguin import load_2D_data
import matplotlib.pyplot as plt
import numpy as np
import geometry as geo
from Interpolate import biInterpolate
from Search import binSearch

def density(x, y, data):
    xVals = data[1]
    yVals = data[2]
    dens = data[3]
    xPoints = binSearch(xVals, 0, len(xVals)-1, x)
    yPoints = binSearch(yVals, 0, len(yVals)-1, y)
    if(xPoints[0] == -np.infty):
        return 0
    elif(xPoints[1] == np.infty):
        return 0
    if(yPoints[0] == -np.infty):
        yPoints = (-1, yPoints[1])
    elif(yPoints[1] == np.infty):
        yPoints = (yPoints[0], -1)
    
    if(xPoints[0] >= len(dens[0]) or xPoints[1] >= len(dens[0])):
        return 0
    if(yPoints[0] >= len(dens) or yPoints[1] >= len(dens)):
        return 0
    lowx = xPoints[0]
    highx = xPoints[1]
    lowy = yPoints[0]
    highy = yPoints[1]
    
    
    squareVals = [dens[lowy,lowx], dens[highy,lowx], dens[lowy,highx], dens[highy,highx]]
    targetCoords = (x,y)
    minCoords = (lowx, lowy)
    maxCoords = (highx, highy)
    den = abs(biInterpolate(targetCoords, squareVals, minCoords, maxCoords))
    return den
    
    
if __name__ == '__main__':
    planetCoords = (1, np.pi)
    planetCoordsCart = geo.sphericalToCartesian(planetCoords, dim = 2)
    print(planetCoordsCart)