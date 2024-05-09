import numpy as np

'''
Algorithm which performs one step in the solution of a system of 1st order ODES to a 5th order approximation, with an adaptive time-step

Inputs:
y0 - Initial conditions (iterable)
xi - Initial position
xf - Final position
maxstep - The largest time step that the algorithm will take, regardless of error calculation
maxerror - Approximate error to be achieved by the solution
fun - function which returns the derivatives of the variables of interest

Returns:
fifth - fifth order solution to the differential equation
error - approximate error in the solution
h - The step size used for this step
hnew - Calculated ideal step size to maintain error budget while also prioritizing speed
'''
def RK45step(y0, xi, xf, maxstep, maxerror, fun):
    h = maxstep
    x = xi
    while(True):
        #Calculating k's
        k1 = h*fun(x, y0)
        k2 = h*fun(x+h/4, y0+k1/4)
        k3 = h*fun(x+3*h/8, y0+3*k1/32+9*k2/32)
        k4 = h*fun(x+12*h/13, y0+1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197)
        k5 = h*fun(x+h, y0 + 439.0*k1/216.0 - 8.0*k2 + 3680*k3/513.0 - 845.0*k4/4104.0)
        k6 = h*fun(x+h/2, y0 - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0)

        #Calculating 5th order approximation and estimated error
        fourth = (25.0*k1/216.0 + 1408*k3/2565.0 + 2197.0*k4/4104.0 - k5/5.0)
        fifth = (16.0*k1/135.0 + 6656.0*k3/12825.0 + 28561.0*k4/56430.0 - 9.0*k5/50.0 + 2.0*k6/55.0)

        error = (fourth - fifth)/h
        for l in range(len(error)):
            if(error[l] < 0):
                error[l] = -error[l]

        #Finding relative error and adapting step size
        estError = 0
        for i in range(len(error)):
            estError += pow(error[i], 2)
        estError = estError**(1/2)

        if(estError > maxerror):
            h = h*0.9*pow(maxerror/estError, 0.2)
            continue

        if(estError <= maxerror):
            hnew = 1.1*h
            if(h > maxstep):
                hnew = maxstep
            break
    return [fifth, error, h, hnew]

'''
Algorithm to return the flow line of a particle in a 2d velocity vector field

Inputs:
xi - x point which the flow line passes through
yi - y point which the flow line passes through
Xs - Collection of grid points in the X direction
Ys - Collection of grid points in the Y direction
maxerror - Desired error for the flow lines
fun - function which returns the x and y velocities at a given point

Returns:
xs - x points which the flow line passes through
ys - y points which the flow line passes through
'''
def flowLine2D(xi, yi, Xs, Ys, maxerror, fun):
    ys = []
    xs = []
    y0 = np.array([xi, yi])
    xf = np.max(Xs)
    yf = np.max(Ys)

    h = abs(xf-xi)
    ys.append(yi)
    xs.append(xi)
    n = 0
    aveCellX = (xf - xi)/len(Xs)
    aveCellY = (yf - yi)/len(Ys)
    while(abs(y0[0]) < xf and abs(y0[1]) < yf and n < 20000):
        correction = RK45step(y0, y0[0], xf, h, maxerror, fun)
        if(abs(correction[0][0]) <= 0.01*aveCellX and abs(correction[0][1]) <= 0.01*aveCellY):
            break
        y0 = y0 + correction[0]
        ys.append(y0[1])
        xs.append(y0[0])
        h = correction[3]
        n = n+1
    y0 = np.array([xi, yi])
    h = abs(xf-xi)
    n = 0
    while(abs(y0[0]) < xf and abs(y0[1]) < yf and n < 20000):
        correction = RK45step(y0, y0[0], xf, h, maxerror, fun)
        if(abs(correction[0][0]) <= 0.01*aveCellX and abs(correction[0][1]) <= 0.01*aveCellY):
            break
        y0 = y0 - correction[0]
        ys.insert(0, y0[1])
        xs.insert(0, y0[0])
        h = correction[3]
        n = n+1
    return (xs, ys)


'''
Algorithm to return the flow line of a particle in a 3d velocity vector field

Inputs:
xi - x point which the flow line passes through
yi - y point which the flow line passes through
zi - z point which the flow line passes through
Xs - Collection of grid points in the X direction
Ys - Collection of grid points in the Y direction
Zs - Collection of grid points in the Z direction
maxerror - Desired error for the flow lines
fun - function which returns the x and y velocities at a given point

Returns:
xs - x points which the flow line passes through
ys - y points which the flow line passes through
zs - z points which the flow line passes through
'''
def flowLine3D(xi, yi, zi, Xs, Ys, Zs, maxerror, fun, maxstep):
    zs = []
    xs = []
    ys = []
    y0 = np.array([xi, yi, zi])
    xf = np.max(np.abs(Xs))
    yf = np.max(np.abs(Ys))
    zf = np.max(np.abs(Zs))
    h = abs(xf-xi)/10000
    zs.append(zi)
    xs.append(xi)
    ys.append(yi)
    n = 0
    aveCellX = abs(xf - xi)/len(Xs)
    aveCellY = abs(yf - yi)/len(Ys)
    aveCellZ = abs(zf - zi)/len(Zs)
    while(abs(y0[0]) <= abs(xf) and abs(y0[1]) <= abs(yf) and abs(y0[2]) <= abs(zf) and n < 20000):
        correction = RK45step(y0, y0[0], xf, h, maxerror, fun)
        if(abs(correction[0][0]) <= 0.001*aveCellX and abs(correction[0][1]) <= 0.001*aveCellY and abs(correction[0][2]) <= 0.001*aveCellZ):
            break
        y0 = y0 + correction[0]
        zs.append(y0[2])
        xs.append(y0[0])
        ys.append(y0[1])
        h = correction[3]
        if(h > maxstep):
            h = maxstep
        n = n+1

    y0 = np.array([xi, yi, zi])
    h = abs(xf-xi)
    n = 0
    while(abs(y0[0]) <= abs(xf) and abs(y0[1]) <= abs(yf) and abs(y0[2]) <= abs(zf) and n < 20000):
        correction = RK45step(y0, y0[0], xf, h, maxerror, fun)
        if(abs(correction[0][0]) <= 0.001*aveCellX and abs(correction[0][1]) <= 0.001*aveCellY and abs(correction[0][2]) <= 0.001*aveCellZ):
            break
        y0 = y0 - correction[0]
        zs.insert(0, y0[2])
        xs.insert(0, y0[0])
        ys.insert(0, y0[1])
        h = correction[3]
        if(h > maxstep):
            h = maxstep
        n = n+1
    return (xs, ys, zs)