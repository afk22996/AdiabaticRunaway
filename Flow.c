#include "Interpolate3DSpherical.c"
#include "Geometry.c"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define xres 288
#define yres 480
#define zres 144
#define coordRes 1000

double *xVals;
double *yVals;
double *zVals;
double *xVel;
double *yVel;
double *zVel;
double planetCoords[3] = {1, M_PI, 0.5*M_PI};
double planetVel[3] = {0, 1, 0};

double* RK5(double* (*func)(double, double*), double* y0, double xi, double h, int size){
    //Allocating memory for calculating k values
    double *k1 = (double *)malloc(size*sizeof(double));
    double *k2 = (double *)malloc(size*sizeof(double));
    double *k3 = (double *)malloc(size*sizeof(double));
    double *k4 = (double *)malloc(size*sizeof(double));
    double *k5 = (double *)malloc(size*sizeof(double));
    double *k6 = (double *)malloc(size*sizeof(double));
    double *tempy = (double *)malloc(size*sizeof(double));
    double *tempk = (double *)malloc(size*sizeof(double));
    double *sol = (double *)malloc(size*sizeof(double));

    tempk = func(xi, y0);
    for(int i = 0; i < size; i++){k1[i] = h*tempk[i];}

    for(int i = 0; i < size; i++){tempy[i] = y0[i] + k1[i]/4.0;}
    tempk = func(xi + h/4.0, tempy);
	for(int i = 0; i < size; i++){k2[i] = h*tempk[i];}

	for(int i = 0; i < size; i++){tempy[i] = y0[i] + 3.0*k1[i]/32.0 + 9.0*k2[i]/32.0;}
	tempk = func(xi + 3.0*h/8.0, tempy);
	for(int i = 0; i < size; i++){k3[i] = h*tempk[i];}

	for(int i = 0; i < size; i++){tempy[i] = y0[i] + 1932.0*k1[i]/2197.0 - 7200.0*k2[i]/2197.0 + 7296.0*k3[i]/2197.0;}
	tempk = func(xi + 12.0*h/13.0, tempy);
	for(int i = 0; i < size; i++){k4[i] = h*tempk[i];}

	for(int i = 0; i < size; i++){tempy[i] = y0[i] + 439.0*k1[i]/216.0 - 8.0*k2[i] + 3680*k3[i]/513.0 - 845*k4[i]/4104.0;}
	tempk = func(xi + h, tempy);
	for(int i = 0; i < size; i++){k5[i] = h*tempk[i];}

	for(int i = 0; i < size; i++){tempy[i] = y0[i] - 8.0*k1[i]/27.0 + 2.0*k2[i] - 3544.0*k3[i]/2565.0 + 1859.0*k4[i]/4104.0 - 11.0*k5[i]/40.0;}
	tempk = func(xi + h/2.0, tempy);
	for(int i = 0; i < size; i++){k6[i] = h*tempk[i];}

	for(int i = 0; i < size; i++){sol[i] = 16.0*k1[i]/135.0 + 6656.0*k3[i]/12825.0 + 28561.0*k4[i]/56430.0 - 9.0*k5[i]/50.0 + 2.0*k6[i]/55.0;}
    return sol;
}

double absmax(double vals[], int len){
	double max = 0;
	for(int i = 0; i < len; i++){
		if(fabs(vals[i]) > max){max = fabs(vals[i]);}
	}
	return max;
}

double absmin(double vals[], int len){
	double min = 1.7e308;
	for(int i = 0; i < len; i++){
		if(vals[i] < 0){vals[i] *= -1;}
		if(vals[i] < min && vals[i] != 0){min = vals[i];}
	}
	return min;
}

double findH(double *r, double* (*func)(double, double*), double *coordX, double *coordY, double *coordZ, int direction){
	double x = r[0];
	double y = r[1];
	double z = r[2];

	double *v = func(x, r);
	double vx = v[0];
	double vy = v[1];
	double vz = v[2];

	int *xPoints = binSearch(coordX, 0, coordRes, x, coordRes);
	int *yPoints = binSearch(coordY, 0, coordRes, y, coordRes);
	int *zPoints = binSearch(coordZ, 0, coordRes, z, coordRes);
	if(xPoints[0] == -1){return 0;}
	else if(xPoints[1] == -1){return 0;}

	if(yPoints[0] == -1){
		yPoints[0] = yres-1;
		yPoints[1] = 0;
	}
	else if(yPoints[1] == -1){
		yPoints[0] = 0;
		yPoints[1] = yres - 1;
	}

	if(zPoints[0] == -1){return 0;}
	if(zPoints[1] == -1){
		if(z > M_PI/2.0){z = M_PI - z;
			zPoints = binSearch(coordZ, 0, coordRes, z, coordRes);
		if(zPoints[0] == -1 || zPoints[1] == -1){return 0;}
		}
	else{
		zPoints[0] = coordRes-1;
		zPoints[1] = coordRes-1;
	}
	}
	int xp = xPoints[direction];
	int yp = yPoints[direction];
	int zp = zPoints[direction];
	double hx, hy, hz;
	if(vx != 0){hx = (coordX[xp] - x)/vx;}
	else{hx = 0;}
	if(vy != 0){hy = (coordY[yp] - y)/vy;}
	else{hy = 0;}
	if(vz != 0){hz = (coordZ[zp] - z)/vz;}
	else{hz = 0;}
	if(hx < 0){hx *= -1;}
	if(hy < 0){hy *= -1;}
	if(hz < 0){hz *= -1;}
	if(hx == 0 && hy == 0 && hz == 0){return 0;}
	double hs[3] = {hx, hy, hz};
	double h = absmin(hs, 3);
	return h;
	}

double myabs(double x){
	if(x < 0){return -x;}
	else{return x;}
}

double **flowLine3D(double xi, double yi, double zi, double Xs[], double Ys[], double Zs[],
	double* (*func)(double, double*), int maxsteps){
	double xf, yf, zf, h, phi, lastPhi, phiInitial;
	double *correction;
	double **sol = (double **)malloc(3*sizeof(double *));
	int n;

	double *xs = (double *)malloc(maxsteps*sizeof(double)+1);
	double *ys = (double *)malloc(maxsteps*sizeof(double)+1);
	double *zs = (double *)malloc(maxsteps*sizeof(double)+1);


	double y0[3] = {xi, yi, zi};

	xf = absmax(Xs, coordRes);
	yf = absmax(Ys, coordRes);
	zf = absmax(Zs, coordRes);

	xs[0] = xi;
	ys[0] = yi;
	zs[0] = zi;

	n = 0;

	phiInitial = cartesianToSpherical(y0, 3)[1];
	lastPhi = phiInitial;
	h = findH(y0, func, Xs, Ys, Zs, 1);
	while((myabs(y0[0]) <= myabs(xf)) && (myabs(y0[1]) <= myabs(yf)) && (myabs(y0[2]) <= myabs(zf))){
		if(h == 0 || h > 5){break;}
		else if(n > maxsteps){break;}
		correction = RK5(func, y0, y0[0], h, 3);
		for(int i = 0; i < 3; i++){
			y0[i] += correction[i];
		}
		xs[n + 1] = y0[0];
		ys[n + 1] = y0[1];
		zs[n + 1] = y0[2];
		h = findH(y0, func, Xs, Ys, Zs, 1);
		phi = cartesianToSpherical(y0, 3)[1];
		if(n > 0){
			if((phi <= phiInitial && lastPhi >= phiInitial) || (phi >= phiInitial && lastPhi <= phiInitial)){break;}
		}
		lastPhi = phi;
		n += 1;
	}
	/*y0[0] = xi;
	y0[1] = yi;
	y0[2] = zi;
	nleft = 0;
	lastPhi = phiInitial;
	h = findH(y0, func, Xs, Ys, Zs, 0);
	while((fabs(y0[0]) <= fabs(xf)) && (fabs(y0[1]) <= fabs(yf)) && (fabs(y0[2]) <= fabs(zf))){
		if(h == 0){break;}
		else if(nleft > maxsteps){break;}
		correction = RK5(func, y0, y0[0], h, 3);
		for(int i = 0; i < 3; i++){
			y0[i] -= correction[i];
		}
		tempx[maxsteps - nleft - 1] = y0[0];
		tempy[maxsteps - nleft - 1] = y0[1];
		tempz[maxsteps - nleft - 1] = y0[2];
		h = findH(y0, func, Xs, Ys, Zs, 0);
		phi = cartesianToSpherical(y0, 3)[1];
		if(nleft > 0){
			if((phi <= phiInitial && lastPhi >= phiInitial) || (phi >= phiInitial && lastPhi <= phiInitial)){break;}
		}
		lastPhi = phi;
		nleft += 1;
	}*/
	sol[0] = xs;
	sol[1] = ys;
	sol[2] = zs;
	printf("%d\n", n);
	return sol;
}

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
	return sphericalToCartesianVelocity(starSphere, v, 3);
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

double* fun(double x, double *y){
	double *sol = (double *)malloc(2*sizeof(double));
	sol[0] = -sin(y[0]);
	sol[1] = y[1];
	return sol;
}

double *starToPlanet(double *sCoords, double *planetCoords){
	double *sCoordsCart = sphericalToCartesian(sCoords, 3);
	double *pCoordsCart = sphericalToCartesian(planetCoords, 3);
	for(int i = 0; i < 3; i++){
		sCoordsCart[i] -= pCoordsCart[i];
	}
	return sCoordsCart;
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
	double guess = 0.01;
	double corot = 0.9969806020039222;
	double sCoords[3] = {corot-guess, 2.0*M_PI/3.0, zVals[zres-2]};
	double *g = starToPlanet(sCoords, planetCoords);
	double **flow = flowLine3D(g[0], g[1], g[2], coordX, coordY, coordZ, vel, 10000);
	for(int i = 0; i < 10001; i++){
		if(flow[0][i] != 0 || flow[1][i] != 0 || flow[2][i] != 0){
			printf("%.15f %.15f %.15f\n", flow[0][i], flow[1][i], flow[2][i]);
		}
	}
	/*double y0[2] = {0.1, 0.2};
	double *s = RK5(fun, y0, 0, 0.1, 2);
	printf("%f, %f\n", s[0], s[1]);*/
}