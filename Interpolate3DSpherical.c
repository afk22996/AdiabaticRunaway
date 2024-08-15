#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int* binSearch(double *vals, int left, int right, double target, int size){
    int middle;
    int *i = (int *)malloc(2*sizeof(int));
    
    while(left < right){
        middle = left + (right - left)/2;
        if(middle >= size-1){break;}
        if(vals[middle] == target){ //Target Value Is In The Array
            i[0] = middle;
            i[1] = middle;
            return i;
        }
        else if(vals[middle+1] == target){
            i[0] = middle + 1;
            i[1] = middle + 1;
            return i;
        }
        else if(vals[middle-1] == target){
            i[0] = middle - 1;
            i[1] = middle - 1;
            return i;
        }
        else if(vals[middle + 1] > target && vals[middle] < target){//Correct Guess With The Next Value Being The Upper Bound
            i[0] = middle;
            i[1] = middle + 1;
            return i;
        }
        else if(vals[middle-1] < target && vals[middle] > target){//Correct Guess With The Previous Value Being The Lower Bound
            i[0] = middle -1;
            i[1] = middle;
            return i;
        }
        else{
            if(vals[middle] < target){left = middle + 1;}//Guess Was Too Low
            else{right = middle - 1;}//Guess Was Too High
        }
    }
    if(right >= size - 1){//Target Was Outside Upper Bound
        i[0] = size - 1;
        i[1] = -1;
    }
    else{//Target Was Outside Lower Bound
        i[0] = -1;
        i[1] = 0;
    }
    return i;
}

double triInterpolate(double *targetCoords, double *cubeVals, double *minCoords, double *maxCoords){
    double c000, c001, c010, c011, c100, c101, c110, c111, x, x0, x1, y, y0, y1, z, z0, z1, xd, yd, zd, c00, c01, c10, c11, c0, c1, c;
    
    //Pulling Out Grid Corners
    c000 = cubeVals[0];
    c001 = cubeVals[1];
    c010 = cubeVals[2];
    c011 = cubeVals[3];
    c100 = cubeVals[4];
    c101 = cubeVals[5];
    c110 = cubeVals[6];
    c111 = cubeVals[7];
    
    //Pullint Target, Lowest, and Highest x,y,z Values
    x = targetCoords[0];
    y = targetCoords[1];
    z = targetCoords[2];
    
    x0 = minCoords[0];
    y0 = minCoords[1];
    z0 = minCoords[2];
    
    x1 = maxCoords[0];
    y1 = maxCoords[1];
    z1 = maxCoords[2];
    
    //Calculating Distance Ratios
    if(x1 != x0){xd = (x-x0)/(x1-x0);}
    else{xd = 0;}
    
    if(y1 != y0){yd = (y-y0)/(y1-y0);}
    else{yd = 0;}
    
    if(z1 != z0){zd = (z-z0)/(z1-z0);}
    else{zd = 0;}
    
    //Linear Interpolation Along x
    c00 = c000*(1-xd) + c100*xd;
    c01 = c001*(1-xd) + c101*xd;
    c10 = c010*(1-xd) + c110*xd;
    c11 = c011*(1-xd) + c111*xd;
    
    //Linear Interpolation Along y
    c0 = c00*(1-yd)+c10*yd;
    c1 = c01*(1-yd)+c11*yd;
    
    //Linear Interpolation Along z
    c = c0*(1-zd) + c1*zd;
    return c;
}
double getData(double *data, int zpoint, int ypoint, int xpoint, int zres, int yres, int xres){
    int idx = 0;
    idx += zpoint*yres*xres;
    idx += ypoint*xres;
    idx += xpoint;
    return data[idx];
}

double interpolate3DSpherical(double *xVals, double *yVals, double *zVals, double *data, double *r, int xres, int yres, int zres){
    double x1, x2, x3;
    int lowx, lowy, lowz, highx, highy, highz;
    
    int *xPoints = (int *)malloc(2*sizeof(int));
    int *yPoints = (int *)malloc(2*sizeof(int));
    int *zPoints = (int *)malloc(2*sizeof(int));
    double *targetCoords = (double *)malloc(3*sizeof(double));
    double *minCoords = (double *)malloc(3*sizeof(double));
    double *maxCoords = (double *)malloc(3*sizeof(double));
    double *cubeVals = (double *)malloc(8*sizeof(double));
    
    x1 = r[0]; //X Position
    x2 = r[1]; //Y Position
    if(x2 < 0){x2 = (x2/(2*M_PI)) + 2*M_PI;}
    x3 = r[2]; //Z Position
    if(x3 > M_PI || x3 < 0){return 0;}
    if(x3 > M_PI/2.0){x3 = M_PI - x3;}
    
    xPoints = binSearch(xVals, 0, xres, x1, xres);
    yPoints = binSearch(yVals, 0, yres, x2, yres);
    zPoints = binSearch(zVals, 0, zres, x3, zres);
    
    if(xPoints[0] == -1){return 0;}
    else if(xPoints[1] == -1){return 0;}
    if(yPoints[0] == -1){
        yPoints[0] = yres-1;
        yPoints[1] = 0;
    }
    
    else if(yPoints[1] == -1){
        yPoints[0] = 0;
        yPoints[1] = yres-1;
    }
    
    if(zPoints[1] == -1){zPoints[1] = zPoints[0];}
    if(zPoints[0] == -1){return 0;}
    lowx = xPoints[0];
    highx = xPoints[1];
    lowy = yPoints[0];
    highy = yPoints[1];
    lowz = zPoints[0];
    highz = zPoints[1];
    
    targetCoords[0] = x1;
    targetCoords[1] = x2;
    targetCoords[2] = x3;
    
    minCoords[0] = xVals[lowx];
    minCoords[1] = yVals[lowy];
    minCoords[2] = zVals[lowz];
    
    maxCoords[0] = xVals[highx];
    maxCoords[1] = yVals[highy];
    maxCoords[2] = zVals[highz];
    
    cubeVals[0] = getData(data, lowz, lowy, lowx, zres, yres, xres);
    cubeVals[1] = getData(data, highz, lowy, lowx, zres, yres, xres);
    cubeVals[2] = getData(data, lowz, highy, lowx, zres, yres, xres);
    cubeVals[3] = getData(data, highz, highy, lowx, zres, yres, xres);
    cubeVals[4] = getData(data, lowz, lowy, highx, zres, yres, xres);
    cubeVals[5] = getData(data, highz, lowy, highx, zres, yres, xres);
    cubeVals[6] = getData(data, lowz, highy, highx, zres, yres, xres);
    cubeVals[7] = getData(data, highz, highy, highx, zres, yres, xres);
    return triInterpolate(targetCoords, cubeVals, minCoords, maxCoords);
    
}
/*int main(void){
    int nPoints = 1000;
    double xVals[nPoints], yVals[nPoints], zVals[nPoints];
    
    for(int i = 0; i < nPoints; i++){
        xVals[i] = 10*i/(nPoints-1);
        yVals[i] = 2*M_PI*i/(nPoints-1);
        zVals[i] = 0.5*M_PI*i/(nPoints-1);
    }
    double *data = (double *)malloc(pow(nPoints, 3)*sizeof(double));
    for(int i = 0; i < nPoints; i++){
        for(int j = 0; j < nPoints; j++){
            for(int k = 0; k < nPoints; k++){
                double r = xVals[i];
                double theta = zVals[k];
                double z = r*cos(theta);
                int idx = k*nPoints*nPoints + j*nPoints + i;
                data[idx] = pow(r, -1.5)*exp(-pow(z, 2.0)/(2.0*pow(0.05,2.0)));
            }
        }
    }
    double r2[3] = {1, M_PI, 0.55*M_PI};
    double intVal2 = interpolate3DSpherical(xVals, yVals, zVals, data, r2, nPoints, nPoints, nPoints);
    printf("%.15f\n", intVal2);}*/