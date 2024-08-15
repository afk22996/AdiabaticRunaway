import imageio.v2 as imageio
from sys import argv

if __name__ == '__main__':
    images = []
    
    nargs = len(argv) - 1
    frames = int(argv[1]) if nargs >= 1 else 101
    frametime = float(argv[2]) if nargs >= 2 else 0.5
    simname = argv[3] if nargs >= 3 else "1J2DIso"
    foldername = argv[4] if nargs >= 4 else "MassFlowVideoFrames"
    plottype = argv[5] if nargs >= 5 else "MassFlow"
    
    readpath = "/home/afkirby/Plots/" + simname + "/" + foldername + "/"
    writepath = "/home/afkirby/Plots/" + simname + "/"
    
    writer = imageio.get_writer(writepath + simname + plottype + ".mp4", fps = 1/frametime)
    for i in range(frames):
        writer.append_data(imageio.imread(readpath + "Orbit" + str(i) + ".png"))
    writer.close()