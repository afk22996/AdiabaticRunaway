import numpy as np


'''Algorithm to convert between spherical coordinates centered on a star to those centered on a planet. Notable assumption is that the planet and the star lie at the same phi coordinate.

Inputs:
starCoords - The coordinates of the object of interest in the star-centric coordinate system
planetCoords - The cooridnates of the planet in the star-centric coordinate system

Output:
Tuple containing r, theta, and phi in the new planet-centric coordinate system

'''
def sphericalStarToPlanet(starCoords, planetCoords):
	rStar = starCoords[0]
	thetaStar = starCoords[1]
	phiStar = starCoords[2]

	rPlanet = planetCoords[0]
	thetaPlanet = planetCoords[1]
	phiPlanet = planetCoords[2]

	rPrime = np.sqrt(np.pow(rStar,2) + np.pow(rPlanet,2) - 2*rStar*rPlanet*np.cos(thetaPlanet - thetaStar))
	thetaPrime = thetaStar-thetaPlanet
	phiPrime = phiStar

	return (rPrime, thetaPrime, phiPrime)


'''Algorithm to convert between spherical and cartesian 3D coordinates

Input:
coords - indexable list of the spherical coordinates (r, theta, phi)

Output - tuple containing cartesian coordinates (x, y, z)'''
def sphericalToCartesian(coords):
	r = coords[0]
	theta = coords[1]
	phi = coords[2]

	x = r*np.sin(phi)*np.cos(theta)
	y = r*np.sin(phi)*np.sin(theta)
	z = r*np.cos(phi)

	return (x, y, z)