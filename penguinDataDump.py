from read_penguin import load_3D_data, cell_center

fpath = "/home/alex/Documents/Coding/Research/Analysis/Data/"
xres = 288
yres = 480
zres = 144
label = "h50_1p10E_e0_PPM4_ave"
frame = 10
data = load_3D_data(fpath, xres, yres, zres, label, frame)

xVals = cell_center(data[1])
yVals = cell_center(data[2])
zVals = cell_center(data[3])
xVels = data[6]
yVels = data[7]
zVels = data[8]

f = open("Xs.txt", 'w')
for i in range(len(xVals)):
	f.write(str(xVals[i]) + "\n")

f = open("Ys.txt", 'w')
for i in range(len(yVals)):
	f.write(str(yVals[i]) + "\n")

f = open("Zs.txt", 'w')
for i in range(len(zVals)):
	f.write(str(zVals[i]) + "\n")

f = open("xVels.txt", 'w')
for k in range(len(zVals)):
	for j in range(len(yVals)):
		for i in range(len(xVals)):
			f.write(str(xVels[k,j,i]) + "\n")

f = open("yVels.txt", 'w')
for k in range(len(zVals)):
	for j in range(len(yVals)):
		for i in range(len(xVals)):
			f.write(str(yVels[k,j,i]) + "\n")

f = open("zVels.txt", 'w')
for k in range(len(zVals)):
	for j in range(len(yVals)):
		for i in range(len(xVals)):
			f.write(str(zVels[k,j,i]) + "\n")