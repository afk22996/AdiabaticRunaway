import numpy as np

def triInterpolate(targetCoords, cubeVals, minCoords, maxCoords):
	c000 = cubeVals[0]
	c001 = cubeVals[1]
	c010 = cubeVals[2]
	c011 = cubeVals[3]
	c100 = cubeVals[4]
	c101 = cubeVals[5]
	c110 = cubeVals[6]
	c111 = cubeVals[7]

	x = targetCoords[0]
	y = targetCoords[1]
	z = targetCoords[2]

	x0 = minCoords[0]
	y0 = minCoords[1]
	z0 = minCoords[2]

	x1 = maxCoords[0]
	y1 = maxCoords[1]
	z1 = maxCoords[2]

	xd = (x-x0)/(x1-x0)
	yd = (y-y0)/(y1-y0)
	zd = (z-z0)/(z1-z0)

	c00 = c000*(1-xd) + c100*xd
	c01 = c001*(1-xd) + c101*xd
	c10 = c010*(1-xd) + c110*xd
	c11 = c011*(1-xd) + c111*xd

	c0 = c00*(1-yd)+c10*yd
	c1 = c01*(1-yd)+c11*yd

	c = c0*(1-zd) + c1*zd
	return c