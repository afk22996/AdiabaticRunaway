#include "Interpolate3DSpherical.c"
#include "Geometry.c"
#include "Flow.c"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#define xres 288
#define yres 480
#define zres 144

double *xVals;
double *yVals;
double *zVals;
double *xVel;
double *yVel;
double *zVel;
double planetCoords[3] = {1, M_PI, 0.5*M_PI};
double planetVel[3] = {0, 1, 0};

double *getVel3D(double *coords, double *planetCoords, double *planetVel, double *xVals, double *yVals, double *zVals, double *xVel, double *yVel, double *zVel){
	double *planetCart = sphericalToCartesian(coords, 3);
	double *starCart = planetCart;
	starCart[0] -= 1;
	double *starSphere = cartesianToSpherical(starCart, 3);
	double vx = interpolate3DSpherical(xVals, yVals, zVals, xVel, starSphere, xres, yres, zres);
	double vy = interpolate3DSpherical(xVals, yVals, zVals, yVel, starSphere, xres, yres, zres) - starSphere[0]*sin(starSphere[2]);
	int sign = (0.5*M_PI >= starSphere[2])?1:-1;
	double vz = sign*interpolate3DSpherical(xVals, yVals, zVals, zVel, starSphere, xres, yres, zres);
	if(coords[2] == 0){vz = 0;}
	double *v = (double *)malloc(3*sizeof(double));
	v[0] = vx;
	v[1] = vy;
	v[2] = vz;
	return v;
}

double *vel(double x, double *y){
	double *coords = cartesianToSpherical(y, 3);
	return getVel3D(coords, planetCoords, planetVel, xVals, yVals, zVals, xVel, yVel, zVel);
}

double *fromFile(char *fname, int size){
	FILE *fp;
	char *line = NULL;
	size_t len = 0;

	fp = fopen(fname, "r");
	double *dat = (double *)malloc(size*sizeof(double));
	if(fp == NULL){
		printf("File not found: %s\n", fname);
		return NULL;
	}
	for(int i = 0; i < size; i++){
		getline(&line, &len, fp);
		dat[i] = atof(line);
	}
	fclose(fp);
	if(line){free(line);}
	return dat;
}

double *linspace(double xi, double xf, int steps){
	double *x = (double *)malloc(steps*sizeof(double));
	for(int i = 0; i < steps; i++){
		x[i] = xi + (xf-xi)*i/(steps-1);
	}
	return x;
}

int main(void){
	xVals = fromFile("Xs.txt", xres);
	yVals = fromFile("Ys.txt", yres);
	zVals = fromFile("Zs.txt", zres);
	xVel = fromFile("xVels.txt", xres*yres*zres);
	yVel = fromFile("yVels.txt", xres*yres*zres);
	zVel = fromFile("zVels.txt", xres*yres*zres);
	double *coordX = linspace(-xVals[xres-1] + 1, xVals[xres-1] + 1, 1000);
	double *coordY = linspace(-xVals[xres-1], xVals[xres-1], 1000);
	double *coordZ = linspace(-cos(zVals[0]), cos(zVals[0]), 1000);
	double **flow;
	flow = flowLine3D(0.01, 0.1, 0, coordX, coordY, coordZ, vel, xres, yres, zres, 1000);
}