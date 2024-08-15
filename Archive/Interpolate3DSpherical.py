import numpy as np

'''Algorithm to interpolate 3D data with a spherical boundary condition assumed which is periodic in Y and mirrored in Z

inputs:
xVals - array of radial cell centers
yVals - array of azimuthal cell centers
zVals - array of polar cell centers
data - data over which to interpolate
r - indexable object which gives the radial, azimuthal, and polar positions at which you want the interpolated data

output:
Interpolated data at the point of interest (returns 0 if outside of the proper z-range)'''
def interpolate3DSpherical(xVals, yVals, zVals, data, r):
    x1 = r[0] #X Position
    x2 = r[1] #Y Position
    if(x2 < 0):
        x2 = x2%(2*np.pi) + 2*np.pi
    x3 = r[2] #Z Position
    if(x3 > np.pi or x3 < 0):
    	return 0
    
    xPoints = binSearch(xVals, 0, len(xVals), x1)
    yPoints = binSearch(yVals, 0, len(yVals), x2)
    zPoints = binSearch(zVals, 0, len(zVals), x3)
    if(xPoints[0] == -np.infty):
        return 0
    elif(xPoints[1] == np.infty):
        return 0
    if(yPoints[0] == -np.infty):
        yPoints = (len(yVals)-1, 0)
    elif(yPoints[1] == np.infty):
        yPoints = (0, len(yVals)-1)
    if(zPoints[1] == np.infty):
        if x3 > np.pi/2:
            x3 = np.pi - x3
            zPoints = binSearch(zVals, 0, len(zVals), x3)
            if(zPoints[0] == -np.infty):
                return 0
            if(zPoints[1] == np.infty):
                return 0
        else:
            zPoints = (-1, -1)
    if(zPoints[0] == -np.infty):
        return 0

    lowx = xPoints[0]
    highx = xPoints[1]
    lowy = yPoints[0]
    highy = yPoints[1]
    lowz = zPoints[0]
    highz = zPoints[1]
    targetCoords = (x1,x2,x3)
    minCoords = (xVals[lowx], yVals[lowy], zVals[lowz])
    maxCoords = (xVals[highx], yVals[highy], zVals[highz])
    cubeVals = [data[lowz,lowy,lowx], data[highz,lowy,lowx], data[lowz,highy,lowx], data[highz,highy,lowx], data[lowz,lowy,highx], data[highz, lowy, highx], data[lowz, highy, highx], data[highz, highy, highx]]
    return triInterpolate(targetCoords, cubeVals, minCoords, maxCoords)


'''Modified binary search which returns a tuple of the two indices between which the target value should lie between. If the value is in the list, the returned tuple contains duplicate values
with both being the index of the target value. If the value is outside of the bounds of the list, a tuple containing either positive or negative infinities is used where the other value is unknown
to commmunicate that uncertainty'''
def binSearch(vals, left, right, target):
    while left < right:
        middle = left + (right - left)//2
        if(middle >= len(vals)-1):
            break

        if(vals[middle] == target): #target value is in array
            return (middle, middle)
        elif(vals[middle + 1] == target):
            return (middle + 1, middle + 1)
        elif(vals[middle - 1] == target):
            return (middle - 1, middle - 1)
        elif(vals[middle +1] > target and vals[middle] < target): #correct guess with the next value being the upper bound
            return (middle, middle+1)

        elif(vals[middle-1] < target and vals[middle] > target): #correct value with the previous value being the lower bound
            return (middle-1, middle)

        else:
            if(vals[middle] < target): #guess was too low
                left = middle + 1

            else: #guess was too high
                right = middle - 1

    if(right > len(vals)-1): #target was outside of upper bound of array
        return(len(vals)-1, np.infty)

    else: #target was outside of lower bound of array
        return(-np.infty,0)

'''
Function which performs a linear interpolation over 3D data given a target coordinate and the values at the grid. This is a 2nd order approximation, and so falls off like d^-2, where d is the distance between lattice points

Inputs:
targetCoords - indexable object which contains the coordinates which you would like the coordinates of
squareVals - indexable object which contains the data values which the target coordinate should fall between. This assumes a cube lattice and should be of the form 
    [minX minY, minZ, minX minY maxZ, minX maxY minZ, minX maxY maxZ, maxX minY minZ, maxX minY maxZ, maxX maxY minZ, maxX maxY maxZ] values of that cube
minCoords - indexable object which contains the lower x/y/z coordinates (i.e. (minCoords[0], minCoords[1], minCoords[2]) would be the coordinates for the bottom left corner out of page corner of the cube lattice)
maxCoords - same as minCoords but contains the higher coordinates

Output:
Approximate data at the target coordinates (2nd order approximation)
'''
def triInterpolate(targetCoords, cubeVals, minCoords, maxCoords):
    #Pulling corner values of the grid
    c000 = cubeVals[0] #left side bottom corner out of page
    c001 = cubeVals[1] #left side top corner out of page
    c010 = cubeVals[2] #left side bottom corner into page
    c011 = cubeVals[3] #left side top corner into page
    c100 = cubeVals[4] #right side bottom corner out of page
    c101 = cubeVals[5] #right side top corner out of page
    c110 = cubeVals[6] #right side bottom corner into page
    c111 = cubeVals[7] #right side top corner into page

    #Pulling target positions
    x = targetCoords[0]
    y = targetCoords[1]
    z = targetCoords[2]

    #Pulling lowest x ,y, and z coordinates
    x0 = minCoords[0]
    y0 = minCoords[1]
    z0 = minCoords[2]

    #Pulling highest x, y, and z coordinates
    x1 = maxCoords[0]
    y1 = maxCoords[1]
    z1 = maxCoords[2]

    #Calculating distance ratios
    if(x1 != x0):
        xd = (x-x0)/(x1-x0)
    else:
        xd = 0
    
    if(y1 != y0):
        yd = (y-y0)/(y1-y0)
    else:
        yd = 0
    if(z1 != z0):
        zd = (z-z0)/(z1-z0)
    else:
        zd = 0

    #Linear interpolation along x
    c00 = c000*(1-xd) + c100*xd
    c01 = c001*(1-xd) + c101*xd
    c10 = c010*(1-xd) + c110*xd
    c11 = c011*(1-xd) + c111*xd

    #Linear interpolation along y
    c0 = c00*(1-yd)+c10*yd
    c1 = c01*(1-yd)+c11*yd

    #Linear interpolation along z
    c = c0*(1-zd) + c1*zd
    return c

if __name__ == '__main__':
    from read_penguin import load_3D_data, cell_center
    from matplotlib import pyplot as plt


    binpath = "C:\\Users\\sidio\\OneDrive\\Desktop\\Coding\\Analysis\\"
    xres = 480
    yres = 864
    zres = 264
    label = "h50_1p10E_e0_PPM4_ave"

    data = load_3D_data(binpath, xres, yres, zres, label, 0)

    rGrid = cell_center(data[1])
    phiGrid = cell_center(data[2])
    thetaGrid = cell_center(data[3])

    thetaPoints = np.linspace(thetaGrid[0], thetaGrid[-1], 1000)
    dens = [interpolate3DSpherical(rGrid, phiGrid, thetaGrid, data[4], (rGrid[240], phiGrid[432], thetaGrid[i])) for i in range(len(thetaGrid))]

    plt.plot(thetaGrid, dens, label = 'interpolated densities', marker = '.', markersize = 5, linewidth = 0)
    #plt.vlines(np.pi/2, 0, 1, linestyles='dashed', label = 'pi/2')
    actualDens = [(rGrid[240]**-1.5)*np.exp((rGrid[240]*np.sin(thetaPoints[i])/0.05**2)*(1/np.sqrt(rGrid[240]**2 + np.cos(thetaPoints[i])**2) - 1/rGrid[240])) for i in range(len(thetaPoints))]
    plt.plot(thetaPoints, np.array(actualDens), label = 'actual densities')
    plt.legend()
    plt.show()