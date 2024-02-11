import numpy as np


'''Algorithm to convert between spherical/polar coordinates centered on a star to those centered on a planet. Notable assumption is that the planet and the star lie at the same phi coordinate.

Inputs:
starCoords - The coordinates of the object of interest in the star-centric coordinate system
planetCoords - The cooridnates of the planet in the star-centric coordinate system
dim - Number of dimensions the conversion is in (2 is polar to cartesian, 3 is spherical)

Output:
Tuple containing r, theta, and phi in the new planet-centric coordinate system

'''
def curvilinearStarToPlanet(starCoords, planetCoords, dim = 3):
	rStar = starCoords[0]
	thetaStar = starCoords[1]

	rPlanet = planetCoords[0]
	thetaPlanet = planetCoords[1]

	rPrime = np.sqrt(np.pow(rStar,2) + np.pow(rPlanet,2) - 2*rStar*rPlanet*np.cos(thetaPlanet - thetaStar))
	thetaPrime = thetaStar-thetaPlanet
	if(dim == 3):
		phiStar = starCoords[2]
		phiPlanet = planetCoords[2]
		phiPrime = phiStar - phiPlanet
		return (rPrime, thetaPrime, phiPrime)
	else:
		return (rPrime, thetaPrime)


'''Algorithm to convert between spherical and cartesian 3D coordinates

Input:
coords - indexable list of the spherical coordinates (r, theta, phi)

Output - tuple containing cartesian coordinates (x, y, z)'''
def sphericalToCartesian(coords, dim = 3):
	r = coords[0]
	theta = coords[1]
	if(dim > 2):
		phi = coords[2]
	else:
		phi = np.pi

	x = r*np.sin(phi)*np.cos(theta)
	y = r*np.sin(phi)*np.sin(theta)
	z = r*np.cos(phi)

	return (x, y, z)

'''Algorithm to transform velocities in a frame of reference centered on a star to the frame of reference of a planet at a known position and with known velocities

Inputs:
starSpeed - Indexable object containing the known velocities in the star-centered reference frame
planetSpeed - Indexable object containing the speed of the planet in the star-centered reference frame
planetPosition - Indexable object containing the position of the planet in the star-centered refrence frame
starPosition - Indexable object containing the position of the point of interest in the star-centered reference frame

Output:
Tuple containing the transformed velocities in 3D'''
def sphericalStarToPlanetVelocity(starSpeed, planetSpeed, planetPosition, starPosition, dim = 3):
	rPrime = sphericalStarToPlanet(starPosition, planetPosition)[0]

	rStar = starPosition[0]
	thetaStar = starPosition[1]

	vRStar = starSpeeds[0]
	vThetaStar = starSpeeds[1]
	

	rPlanet = planetPosition[0]
	thetaPlanet = planetPosition[1]

	vRPlanet = planetSpeed[0]
	vThetaPlanet = planetSpeed[1]
	

	if(dim > 2):
		vPhiStar = starSpeeds[2]
		vPhiPlanet = planetSpeed[2]

	else:
		vPhiStar = 0
		vPhiPlanet = 0


	vRPrime = (rStar - rPlanet*np.cos(thetaPlanet - thetaStar))*vRStar/rPrime + rStar*rPlanet*np.sin(thetaPlanet - thetaStar)*(vThetaPlanet - vThetaStar)/rPrime
	vThetaPrime = vThetaPlanet - vThetaStar
	vPhiPrime = vPhiPlanet - vPhiStar

	return (vRPrime, vThetaPrime, vPhiPrime)