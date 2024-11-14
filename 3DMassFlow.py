from getQuantity import getVel3D, getDen3D
from read_penguin import load_3D_data, cell_center
import Geometry as geo
import numpy as np
from Numerics import doubleSimpson
from matplotlib import pyplot as plt
planetCoords = (1,np.pi,np.pi/2)
planetVel = (0,1,0)

def getvr3D(coords, data):
    vel = np.array(getVel3D(coords, planetCoords, planetVel, data, corot = True, cart = True))
    vr = -vel[0]*np.cos(coords[1])*np.sin(coords[2]) - vel[1]*np.sin(coords[1])*np.sin(coords[2]) + vel[2]*np.cos(coords[2])
    return vr

def getMDot3DPlus(r, nt, nphi, data):
    def func(phi, theta):
        coords = np.array([r, phi, theta])
        vr = getvr3D(coords, data)
        den = getDen3D(r, phi, theta, data, planetCoords)
        mdot =  -vr*den*(r**2)*np.sin(theta)
        if(mdot > 0):
            return mdot
        else:
            return 0
    return 2*doubleSimpson(func, (0, 2*np.pi), (0, np.pi/2), nt, nphi)

xres3d = 288
yres3d = 576
zres3d = 168
q = 0.00003
h = 0.05
qth = q/h**3
rb3 = qth*h
rh3 = (1/3)**(1/3)*qth**(-2/3)*rb3
choksi = 3.5*qth**2*h**3
data = load_3D_data("/scratch/afkirby/", xres3d, yres3d, zres3d, "h50_1p10E_e0_PPM4_ave", 13)
'''observed = getMDot3DPlus(min(rb3, rh3), 101, 101, data)
bondi = rb3**2*h

print(observed/bondi, choksi/bondi)
print(abs(observed - choksi)/abs(choksi))'''
N = 300
vr = np.zeros((N, N))
phis = np.linspace(0, 2*np.pi, N)
thetas = np.linspace(0, np.pi/2, N)

for i in range(N):
    for j in range(N):
        coords = (rb3, phis[i], thetas[j])
        vr[i,j] = getvr3D(coords, data)

plt.pcolormesh(phis, thetas, vr, cmap = 'hsv')
plt.colorbar()
plt.xlabel("Longitude (rad)")
plt.ylabel("Latitude (rad)")
plt.title("Radial Velocity")
plt.savefig("vrtestFung.png")