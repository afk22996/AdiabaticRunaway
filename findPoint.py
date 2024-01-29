from Interpolate import triInterpolate
from Search import binSearch
import numpy as np

def findPoint(targetCoords, data, geom="sphere", dim = 3):
	if dim == 3:
		rdata = np.ndarray(9)
		rdata[0] = data[0] #Simulation Frame
		rdata[1] = targetCoords[0] #X Position
		rdata[2] = targetCoords[1] #Y Position
		rdata[3] = targetCoords[2] #Z Position


		#Finding nearest points in x,y, and z directions
		xpoints = binSearch(data[1], 0, len(data[1]), rdata[1])
		if(xpoints[0] == -np.pi):
			xpoints = (xpoints[1], xpoints[1])

		if(xpoints[1] == np.pi):
			xpoints = (xpoints[0], xpoints[0])

		ypoints = binSearch(data[2], 0, len(data[2]), rdata[2])
		if(geom == "sphere"): #Applying spherical boundary conditions
			if (ypoints[0] == -np.infty):
				ypoints = (len(data[2])-1, ypoints[1])

			if (ypoints[1] == np.infty):
				ypoints = (ypoints[0], 0)

		ztarget = rdata[3]
		if(rdata[3] < 0): #Applying mirror z symmetry
			ztarget = -ztarget
		zpoints = binSearch(data[3], 0, len(data[3]), ztarget)
		if(geom == "sphere"): #Applying spherical boundary conditions
			if (zpoints[0] == -np.infty):
				zpoints = (len(data[3])-1, zpoints[1])

			if (zpoints[1] == np.infty):
				zpoints = (zpoints[0], 0)

		#Making Cube Lattice
		cubeVals = np.ndarray(8)
		lowx = xpoints[0]
		highx = xpoints[1]
		lowy = ypoints[0]
		highy = ypoints[1]
		lowz = zpoints[0]
		highz = zpoints[1]

		minCoords = (lowx, lowy, lowz)
		maxCoords = (highx, highy, highz)

		#Interpolating Density
		cubeVals[0] = data[4][lowz, lowy, lowx]
		cubeVals[1] = data[4][highz, lowy, lowx]
		cubeVals[2] = data[4][lowz, highy, lowx]
		cubeVals[3] = data[4][highz, highy, lowx]
		cubeVals[4] = data[4][lowz, lowy, highx]
		cubeVals[5] = data[4][highz, lowy, highx]
		cubeVals[6] = data[4][lowz, highy, highx]
		cubeVals[7] = data[4][highz, highy, highx]

		rdata[4] = triInterpolate(targetCoords, cubeVals, mincoords, maxcoords)

		#Interpolating Pressure
		cubeVals[0] = data[5][lowz, lowy, lowx]
		cubeVals[1] = data[5][highz, lowy, lowx]
		cubeVals[2] = data[5][lowz, highy, lowx]
		cubeVals[3] = data[5][highz, highy, lowx]
		cubeVals[4] = data[5][lowz, lowy, highx]
		cubeVals[5] = data[5][highz, lowy, highx]
		cubeVals[6] = data[5][lowz, highy, highx]
		cubeVals[7] = data[5][highz, highy, highx]

		rdata[5] = triInterpolate(targetCoords, cubeVals, mincoords, maxcoords)

		#Interpolating X-Velocity
		cubeVals[0] = data[6][lowz, lowy, lowx]
		cubeVals[1] = data[6][highz, lowy, lowx]
		cubeVals[2] = data[6][lowz, highy, lowx]
		cubeVals[3] = data[6][highz, highy, lowx]
		cubeVals[4] = data[6][lowz, lowy, highx]
		cubeVals[5] = data[6][highz, lowy, highx]
		cubeVals[6] = data[6][lowz, highy, highx]
		cubeVals[7] = data[6][highz, highy, highx]

		rdata[6] = triInterpolate(targetCoords, cubeVals, mincoords, maxcoords)

		#Interpolating Y-Velocity
		cubeVals[0] = data[7][lowz, lowy, lowx]
		cubeVals[1] = data[7][highz, lowy, lowx]
		cubeVals[2] = data[7][lowz, highy, lowx]
		cubeVals[3] = data[7][highz, highy, lowx]
		cubeVals[4] = data[7][lowz, lowy, highx]
		cubeVals[5] = data[7][highz, lowy, highx]
		cubeVals[6] = data[7][lowz, highy, highx]
		cubeVals[7] = data[7][highz, highy, highx]

		rdata[7] = triInterpolate(targetCoords, cubeVals, mincoords, maxcoords)

		#Interpolating Z-Velocity
		cubeVals[0] = data[8][lowz, lowy, lowx]
		cubeVals[1] = data[8][highz, lowy, lowx]
		cubeVals[2] = data[8][lowz, highy, lowx]
		cubeVals[3] = data[8][highz, highy, lowx]
		cubeVals[4] = data[8][lowz, lowy, highx]
		cubeVals[5] = data[8][highz, lowy, highx]
		cubeVals[6] = data[8][lowz, highy, highx]
		cubeVals[7] = data[8][highz, highy, highx]

		rdata[8] = triInterpolate(targetCoords, cubeVals, mincoords, maxcoords)

		return rdata
	return