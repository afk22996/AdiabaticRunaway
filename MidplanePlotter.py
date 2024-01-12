from read_penguin import load_3D_data
import numpy as np
import matplotlib.pyplot as plt

#Binary File Paramters, the binary file should be in the directory corresponding to path with the title "binary_XresxYresxZres_label_time"
Xres = 384
Yres = 768
Zres = 216
path = "/scratch/afkirby/3DCPDAdiabatic/"
label = "h50_1p10E_e0_PPM4"
time = 20

data = load_3D_data(path, Xres, Yres, Zres, label, time)
xpos = data[1] #R positions
ypos = data[2] #Theta positions
zpos = data[3] #Phi positions
density = data[4] #3d array of densities indicies are in reverse order, so density[1,2,3] is the density at z = 1, y = 2, and x = 3

filename = "AdiabaticMidplaneFinal.png"
title = "Final Midplane Density"
xlabel = "Radius"
ylabel = "Azimutal Angle (rad)"

plt.pcolor(xpos, ypos, density[0]) #taking only the densities at z = 0 (midplane)
plt.title(title)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.colorbar()


plt.savefig(filename)