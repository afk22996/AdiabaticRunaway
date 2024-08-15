import numpy as np
import read_penguin as rd
import matplotlib.pyplot as plt

def get_vr(vr_star,vp_star,r_star,phi_star,r_planet,phi_planet): #get vr in planet frame from vr in star frame
    cos_star = np.cos(phi_star)
    sin_star = np.sin(phi_star)

    vr_planet = (vr_star*cos_star - vp_star*sin_star)*np.cos(phi_planet) + (vr_star*sin_star + vp_star*cos_star)*np.sin(phi_planet)
    return vr_planet

def planet_to_star(rad,azi,planet_r,planet_phi): #convert from planet-centric polar to star-centric polar
    x_star = rad*np.cos(azi)+planet_r*np.cos(planet_phi)
    y_star = rad*np.sin(azi)+planet_r*np.sin(planet_phi)

    r_star = np.sqrt(x_star*x_star + y_star*y_star)
    phi_star = np.arctan2(y_star,x_star)
    return r_star, phi_star

def Simpson_quadrature(func,h): #just the area under a curve
    return (h/3.0)*(np.sum(func)+np.sum(func[1:-1])+2.0*np.sum(func[1::2]))

def get_M_dot(dat,planet_r,planet_phi,rad,N,debug=False):
    # inputs are:
    # dat directly from penguin
    # planet_r the radial position of the planet
    # planet_phi the azimuthal positon of the planet
    # rad the radius of the circule that we are drawing around the planet
    # N the number of sampling points (must be odd for Simpson quadrature)

    azi = np.linspace(0.0,2.0*np.pi,num=N)
    r_star, phi_star = planet_to_star(rad,azi,planet_r,planet_phi)

    r_grid = rd.cell_center(dat[1]) #don't forget to create the grid from cell boundaries
    phi_grid = rd.cell_center(dat[2])

    vr_star = np.zeros(N)
    vp_star = np.zeros(N)
    sigma = np.zeros(N)
    for i in range(N):
        sigma[i]   = rd.interpolate_2D_periodicY((r_star[i],phi_star[i]),r_grid,phi_grid,dat[3],dat[2][0],dat[2][-1])
        vr_star[i] = rd.interpolate_2D_periodicY((r_star[i],phi_star[i]),r_grid,phi_grid,dat[5],dat[2][0],dat[2][-1])
        vp_star[i] = rd.interpolate_2D_periodicY((r_star[i],phi_star[i]),r_grid,phi_grid,dat[6],dat[2][0],dat[2][-1])

    vr_planet = get_vr(vr_star,vp_star-r_star,r_star,phi_star,rad,azi)
    func = rad * sigma * vr_planet

    if debug:
        for i in range(N):
            print(i, r_star[i],phi_star[i], sigma[i],vr_star[i],vp_star[i],vr_planet[i])

    return Simpson_quadrature(func,azi[1]-azi[0])

def plot_M_dot_vs_rad(dat,planet_r,planet_phi, frame, ylim):
    # inputs are:
    # dat directly from penguin
    # planet_r the radial position of the planet
    # planet_phi the azimuthal positon of the planet

    rad = np.linspace(0.0,0.3,100)
    M_dot = np.zeros(100)

    for i in range(100):
        M_dot[i]=get_M_dot(dat,planet_r,planet_phi,rad[i],i*10+11)

    plt.plot(rad,M_dot)
    plt.title("Orbit " + str(frame))
    plt.ylim(ylim)
    plt.xlabel("Radius (planet a)")
    plt.show()
    return

def return_M_dot_rad(dat, planet_r, planet_phi):
    rad = np.linspace(0.0,0.3,100)
    M_dot = np.zeros(100)

    for i in range(100):
        M_dot[i]=get_M_dot(dat,planet_r,planet_phi,rad[i],i*10+11)
    return (rad, M_dot)
    