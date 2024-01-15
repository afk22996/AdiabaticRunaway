import numpy as np

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
	xd = (x-x0)/(x1-x0)
	yd = (y-y0)/(y1-y0)
	zd = (z-z0)/(z1-z0)

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
	def f(x,y,z):
		return x**5 + y**5 + z**5

	targetCoords = (1.5, 1.5, 1.5)

	cubeVals = np.ndarray(8)
	low = 1.4999
	high = 1.5001
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
	print("Interpolation: " + str(interpolate))
	print("Actual Value: " + str(actual))
	print("Relative Error: " + str(abs(actual-interpolate)/actual))