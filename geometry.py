import numpy as np

'''
Function to transform cartesian coordinates to spherical/polar coordinates

Inputs:
coords - indexable array representing the cartesian coordinates you want to transform
dim - dimensions of the transformation (2 = polar, 3 = spherical)

Output:
Tuple containing the radius, polar angle, and azimuthal angle. If dim = 2, phi is set to 0
'''
def cartesianToSpherical(coords, dim = 2):
    x = coords[0]
    y = coords[1]
    if(dim > 2):
        z = coords[2]
    else:
        z = 0
    r = np.sqrt(x**2 + y**2 + z**2)
    if(y == 0):
        if x > 0:
            theta = 0
        else:
            theta = -np.pi
    else:
        theta = np.arctan2(y,x)
    phi = np.arccos(z/r)
    if(theta < 0):
        theta = 2*np.pi + theta
    if(dim == 2):
        return (r,theta)
    return(r, theta, phi)

'''
Function to transform spherical/polar coordinates to cartesian coordinates

Inputs:
coords - indexable array representing the polar/spherical coordinates you want to transform
dim - dimensions of the transformation (2 = polar, 3 = spherical)

Output:
Tuple containing the x, y, and z coordinates at the given polar/spherical point. If dim = 2, z is set to 0
'''
def sphericalToCartesian(coords, dim = 2):
	r = coords[0]
	theta = coords[1]
	if(dim > 2):
		phi = coords[2]
	else:
		phi = np.pi/2
	x = r*np.sin(phi)*np.cos(theta)
	y = r*np.sin(phi)*np.sin(theta)
	z = r*np.cos(phi)
	if(dim == 2):
	    return (x,y)
	return (x, y, z)

'''
Function to transform velocities in a spherical/polar basis to a cartesian one. Unlike the coordinate transformations above, this changes the magnitude of each component of the vector, but leaves the overall vector magnitude unchanged.

Inputs:
coords - indexable array representing the polar/spherical coordinates you want to transform the velocities at
velocities - the polar/spherical velocities which you would like to transform
dim - dimensions of the transformation (2 = polar, 3 = spherical).

Output:
Tuple containing the x and y velocities
'''
def sphericalToCartesianVelocity(coords, velocities, dim = 2):
    r = coords[0]
    vr = velocities[0]
    phi = coords[1]
    vphi = velocities[1]
    if(dim > 2):
        theta = coords[2]
        vtheta = velocities[2]
    else:
        theta = np.pi/2
        vtheta = 0

    vx = vr*np.sin(theta)*np.cos(phi) + r*vtheta*np.cos(theta)*np.cos(phi) - r*vphi*np.sin(theta)*np.sin(phi)
    vy = vr*np.sin(theta)*np.sin(phi) + r*vtheta*np.cos(theta)*np.sin(phi) + r*vphi*np.sin(theta)*np.cos(phi)
    vz = vr*np.cos(theta) - r*theta*np.sin(theta)
    if(dim > 2):
        return(vx, vy, vz)

    return (vx, vy)

'''
Function to transform velocities in a  cartesian basis to a spherical/polar one. Unlike the coordinate transformations above, this changes the magnitude of each component of the vector, but leaves the overall vector magnitude unchanged.
This also currently only works for polar transformations, but will be changed to work in 3D as well.

Inputs:
coords - indexable array representing the cartesian coordinates you want to transform the velocities at
velocities - the cartesian velocities which you would like to transform
dim - dimensions of the transformation (2 = polar, 3 = spherical). Currently unused but will be added soon

Output:
Tuple containing the radial and angular velocities
'''
def cartesianToSphericalVelocity(coords, velocities, dim = 2):
    x = coords[0]
    y = coords[1]
    vx = velocities[0]
    vy = velocities[1]
    if(dim > 2):
        z = coords[2]
        vz = velocities[2]
    else:
        z = 0
        vz = 0
    r = np.sqrt(x**2 + y**2 + z**2)
    if(r == 0):
        return (0,0,0)
    vr = (x*vx + y*vy + z*vz)/r
    vphi = (x*vy - y*vx)/(x**2 + y**2)
    vtheta = (z*(x*vx + y*vy) - (x**2 + y**2)*vz)/(r**2*np.sqrt(x**2 + y**2))
    if(dim > 2):
        return (vr, vphi, vtheta)
    return (vr, vphi)


#Testing
if __name__ == '__main__':
	from random import random
	from matplotlib import pyplot as plt
	differences = []
	for i in range(10000):
	    initial = (10*random(), 2*np.pi*random())
	    initialv = (10*random(), 2*np.pi*random())
	    polar = cartesianToSpherical(initial, dim = 2)
	    polarv = cartesianToSphericalVelocity(initial, initialv, dim=2)
	    finalv = sphericalToCartesianVelocity(polar, polarv, dim=2)
	    differences.append(np.sqrt((finalv[0] - initialv[0])**2 + (finalv[1] - initialv[1])**2))
	plt.plot(differences)
	plt.show()