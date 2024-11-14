import numpy as np

def RK5Step(y0, xi, xf, maxstep, fun):
    h = maxstep
    x = xi
    #Calculating k's
    k1 = h*fun(x, y0)
    k2 = h*fun(x+h/4, y0+k1/4)
    k3 = h*fun(x+3*h/8, y0+3*k1/32+9*k2/32)
    k4 = h*fun(x+12*h/13, y0+1932*k1/2197 - 7200*k2/2197 + 7296*k3/2197)
    k5 = h*fun(x+h, y0 + 439.0*k1/216.0 - 8.0*k2 + 3680*k3/513.0 - 845.0*k4/4104.0)
    k6 = h*fun(x+h/2, y0 - 8.0*k1/27.0 + 2.0*k2 - 3544.0*k3/2565.0 + 1859.0*k4/4104.0 - 11.0*k5/40.0)
    
    fifth = (16.0*k1/135.0 + 6656.0*k3/12825.0 + 28561.0*k4/56430.0 - 9.0*k5/50.0 + 2.0*k6/55.0)

    return np.array([fifth, h])