import numpy as np
from matplotlib import pyplot as plt

'''
Function which performs a linear interpolation over 2D data given a target coordinate and the values at the grid. This is a 2nd order approximation, and so falls off like d^-2, where d is the distance between lattice points

Inputs:
targetCoords - indexable object which contains the coordinates which you would like the coordinates of
squareVals - indexable object which contains the data values which the target coordinate should fall between. This assumes a square lattice and should be of the form [bottom left, top left, bottom right, top right] values of that square
minCoords - indexable object which contains the lower x/y coordinates (i.e. (minCoords[0], minCoords[1]) would be the coordinates for the bottom left corner of the square lattice)
maxCoords - same as minCoords but contains the higher coordinates

Output:
Approximate data at the target coordinates (2nd order approximation)
'''
def biInterpolate(targetCoords, squareVals, minCoords, maxCoords):
    c00 = squareVals[0] #Bottom left corner of square
    c01 = squareVals[1] #Top left corner of square
    c10 = squareVals[2] #Bottom right corner of square
    c11 = squareVals[3] #Top right corner of square
    
    #Pulling targt positions
    x = targetCoords[0]
    y = targetCoords[1]
    
    #Pulling lowest x, and y coordinates
    x0 = minCoords[0]
    y0 = minCoords[1]
    
    #Pulling highest x, and y coordinates
    x1 = maxCoords[0]
    y1 = maxCoords[1]
    
    #Linear Interpolation in X
    if(x1 != x0):
        c0 = (x1-x)/(x1-x0)*c00 + (x-x0)/(x1-x0)*c10
        c1 = (x1-x)/(x1-x0)*c01 + (x-x0)/(x1-x0)*c11
    else:
        c0 = c00
        c1 = c01
    
    #Linear Interpolation in Y
    if(y1 != y0):
        c = (y1-y)/(y1-y0)*c0 + (y-y0)/(y1-y0)*c1
    else:
        c = c0
    
    return c
    
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

#Testing
if __name__ == '__main__':
	def f(x,y,z):
		return x**2 + y**2 + z**2

	targetCoords = (1.5, 1.5, 1.5)
	x = []
	y = []
	for i in range(0, 1000000):
		cubeVals = np.ndarray(8)
		low = 1.5 - (i+1)/1000000
		high = 1.5 + (i+1)/1000000
		cubeVals[0] = f(low,low,low)
		cubeVals[1] = f(low,low,high)
		cubeVals[2] = f(low,high,low)
		cubeVals[3] = f(low,high,high)
		cubeVals[4] = f(high,low,low)
		cubeVals[5] = f(high,low,high)
		cubeVals[6] = f(high,high,low)
		cubeVals[7] = f(high,high,high)

		minCoords = (1,1,1)
		maxCoords = (2,2,2)

		interpolate = triInterpolate(targetCoords, cubeVals, minCoords, maxCoords)
		actual = f(1.5,1.5,1.5)
		relerr = abs(actual-interpolate)/actual
		x.append(high-low)
		y.append(relerr)
	plt.plot(x, y)
	plt.xlabel("Grid Size")
	plt.ylabel("Relative Error")
	plt.title("Trilinear Interpolation Error")
	plt.xscale("log")
	plt.yscale("log")
	plt.show()
	'''print("Interpolation: " + str(interpolate))
	print("Actual Value: " + str(actual))
	print("Relative Error: " + str(abs(actual-interpolate)/actual))'''