from matplotlib import pyplot as plt
import numpy as np

datafile = open("/home/alex/Documents/Coding/Research/AdiabaticRunaway/test.txt", 'r')
zs = []
xs = []
ys = []
lines = datafile.readlines()
for line in lines:
    line = line.split(" ")
    xs.append(float(line[0]))
    ys.append(float(line[1]))
    zs.append(float(line[2]))

ax = plt.figure().add_subplot(projection = '3d')
ax.plot(xs, ys, zs)
ax.set_zlim(0, 0.1)
plt.show()