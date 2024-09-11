#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double *cartesianToSpherical(double *coords, int dim){
	double x,y,z,r,theta,phi;
	double *sCoords = (double *)malloc(dim*sizeof(double));

	x = coords[0];
	y = coords[1];
	if(dim > 2){z = coords[2];}
	else{z = 0;}

	r = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
	if(y == 0){
		if(x > 0){theta = 0;}
		else{theta = -M_PI;}
	}
	else{theta = atan2(y,x);}
	phi = acos(z/r);
	if(theta < 0){theta = 2*M_PI + theta;}
	if(r <= 1e-10){
		r = 0;
		theta = 0;
		phi = 0;
	}
	sCoords[0] = r;
	sCoords[1] = theta;
	if(dim > 2){
		sCoords[2] = phi;
	}
	return sCoords;
}

double *sphericalToCartesian(double *coords, int dim){
	double r, theta, phi, x, y, z;
	double *cCoords = (double *)malloc(dim*sizeof(double));

	r = coords[0];
	phi = coords[1];
	if(dim > 2){theta = coords[2];}
	else{theta = M_PI/2.0;}
	x = r*sin(theta)*cos(phi);
	y = r*sin(theta)*sin(phi);
	z = r*cos(theta);
	cCoords[0] = x;
	cCoords[1] = y;
	if(dim > 2){cCoords[2] = z;}
	return cCoords;
}

double *sphericalToCartesianVelocity(double *coords, double *velocities, int dim){
	double r,vr,phi,vphi,theta,vtheta,vx,vy,vz;
	double *cVels = (double *)malloc(dim*sizeof(double));

	r = coords[0];
	vr = velocities[0];
	phi = coords[1];
	vphi = velocities[1];
	if(dim > 2){
		theta = coords[2];
		vtheta = velocities[2];
	}
	else{
		theta = M_PI/2.0;
		vtheta = 0;
	}

	vx = vr*sin(theta)*cos(phi) + vtheta*cos(theta)*cos(phi) - vphi*sin(phi);
	vy = vr*sin(theta)*sin(phi) + vtheta*cos(theta)*sin(phi) + vphi*cos(phi);
	vz = vr*cos(theta) - vtheta*sin(theta);

	cVels[0] = vx;
	cVels[1] = vy;
	if(dim > 2){cVels[2] = vz;}

	return cVels;
}

double *cartesianToSphericalVelocity(double *coords, double *velocities, int dim){
	double x,y,z,vx,vy,vz,r,vr,vphi,vtheta;
	double *sVels = (double *)malloc(dim*sizeof(double));

	x = coords[0];
	y = coords[1];
	vx = velocities[0];
	vy = velocities[1];
	if(dim > 2){
		z = coords[2];
		vz = velocities[2];
	}
	else{
		z = 0;
		vz = 0;
	}
	r = sqrt(pow(x,2) + pow(y,2) + pow(z,2));
	if(r == 0){
		sVels[0] = 0;
		sVels[1] = 0;
		sVels[2] = 0;
		return sVels;
	}
	vr = (x*vx + y*vy + z*vz)/r;
	vphi = ((x*vy - y*vx)/(pow(x,2) + pow(y,2)))*r;
	vtheta = ((z*(x*vx + y*vy) - vz*(pow(x,2) + pow(y,2)))/(pow(r,2)*sqrt(pow(x,2) + pow(y,2))))*r*sin(z/r);

	sVels[0] = vr;
	sVels[1] = vphi;
	if(dim > 2){sVels[2] = vtheta;}
	return sVels;
}

int main(void){
	double coords[3] = {0.1, 0.2, 0.3};
	double vels[3] = {-1, 0.32, 2.14};
	double *tCoords = cartesianToSphericalVelocity(coords, vels, 3);
	printf("%.15f, %.15f, %.15f\n", tCoords[0], tCoords[1], tCoords[2]);
	return 0;
}