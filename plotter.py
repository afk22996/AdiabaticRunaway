from matplotlib import pyplot as plt
import numpy as np


plt.figure()
datafile = open("/home/afkirby/Corotation10.txt", 'r')
zs = []
xs = []
lowzs = []
lines = datafile.readlines()
for line in lines:
    line = line.split(" ")
    zs.append(float(line[0]))
    lowzs.append(float(line[-1]))
    xs.append(0.5*(float(line[1]) + float(line[2])))
plt.plot(zs, lowzs)
plt.ylabel("Lowest Z")
plt.xlabel("Z")
plt.savefig("/home/afkirby/lowestz10.png")

plt.figure()
plt.plot(zs, xs)
plt.ylabel("Horseshoe Half-Widths")
plt.xlabel("Z")
plt.savefig("/home/afkirby/halfwidths10.png")
