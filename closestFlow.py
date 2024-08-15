import Geometry as geo
from getQuantity import getVel3D
from Flow import flowLine3D, isHorseshoe
import numpy as np
from read_penguin import load_3D_data, cell_center
from Numerics import finDiff, NR

planetCoords = (1, np.pi, np.pi/2)
planetVel = (0, 1, 0)

xres = 288
yres = 480
zres = 144

orbit = 10
data = load_3D_data("/scratch/afkirby/3DRsVariation/rs001/", xres, yres, zres, "h50_1p10E_e0_PPM4_ave", orbit)


def vel(x,y):
    coords = geo.cartesianToSpherical(y, dim = 3)
    v = getVel3D(coords, planetCoords, planetVel, data, cart = True)
    return v

def starToPlanet(sCoords, planetCoords):
    sCoordsCart = geo.sphericalToCartesian(sCoords, dim = 3)
    pCoordsCart = np.array(sCoordsCart) - np.array(geo.sphericalToCartesian(planetCoords, dim = 3))
    return pCoordsCart


if __name__ == "__main__":
    f = open("/home/afkirby/Corotation" + str(orbit) + "Adi.txt", 'w+')
    
    rGrid = cell_center(data[1])
    thetaGrid = cell_center(data[3])
    
    xp,yp,zp = geo.sphericalToCartesian(planetCoords, dim = 3)
    
    coordX = np.linspace(-rGrid[-1] - xp, rGrid[-1] - xp, 1000)
    coordY = np.linspace(-rGrid[-1] - yp, rGrid[-1] - yp, 1000)
    coordZ = np.linspace(-np.cos(thetaGrid[0]) - zp, np.cos(thetaGrid[0]) - zp, 1000)
    
    gStep = abs(rGrid[xres//2] - rGrid[xres//2 - 1])
    error = 1e-8
    flowStep = 0.01
    for i in range(1, len(thetaGrid) + 1):
        guess = gStep
        lastGuess = guess
        
        z = np.cos(thetaGrid[-i])
        def velMag(r):
            phi = 2*np.pi/3
            theta = thetaGrid[-i]
            coords = starToPlanet((r, phi, theta), planetCoords)
            v = vel(0, coords)
            return np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
        
        corot = NR(velMag, 1, 1e-3, 1e-10)
        if(corot < rGrid[0] or corot > rGrid[1]):
            corot = 1
            print(str(z) + ": corotation region was not found")
        
        horseshoe = False
        whileFlag = True
        upoints = []
        lpoints = []
        
        while(whileFlag):
            guCoords = (corot-guess, 2*np.pi/3, thetaGrid[-i])
            gux, guy, guz = starToPlanet(guCoords, planetCoords)
            uflow = flowLine3D(gux, guy, guz, coordX, coordY, coordZ, error, vel, flowStep)
            
            uhorseshoe = isHorseshoe(uflow, planetCoords)
            if(uhorseshoe):
                upoints.append(guess)
            
            
            glCoords = (corot+guess, 4*np.pi/3, thetaGrid[-i])
            glx, gly, glz = starToPlanet(glCoords, planetCoords)
            lflow = flowLine3D(glx, gly, glz, coordX, coordY, coordZ, error, vel, flowStep)
            
            lhorseshoe = isHorseshoe(lflow, planetCoords)
            if(lhorseshoe):
                lpoints.append(guess)
            guess = guess + gStep
            if(guess < 0 or guess >= planetCoords[0]):
                whileFlag = False
        
        s = str(z)
        if(len(upoints) > 0):
            gufCoords = (corot-np.max(upoints), 2*np.pi/3, thetaGrid[-i])
            gxuf, gyuf, gzuf = starToPlanet(gufCoords, planetCoords)
            largestuFlow = flowLine3D(gxuf, gyuf, gzuf, coordX, coordY, coordZ, error, vel, flowStep)

            giuCoords = (corot-np.min(upoints), 2*np.pi/3, thetaGrid[-i])
            gxui, gyui, gzui = starToPlanet(giuCoords, planetCoords)
            smallestuFlow = flowLine3D(gxui, gyui, gzui, coordX, coordY, coordZ, error, vel, flowStep)

            s = s + (" " + str(np.min(upoints)) + " " + str(np.min(smallestuFlow[2])) + " " + str(np.max(upoints)) + " " + str(np.min(largestuFlow[2])))
        if(len(lpoints) > 0):
            glfCoords = (corot-np.max(lpoints), 2*np.pi/3, thetaGrid[-i])
            gxlf, gylf, gzlf = starToPlanet(glfCoords, planetCoords)
            largestlFlow = flowLine3D(gxlf, gylf, gzlf, coordX, coordY, coordZ, error, vel, flowStep)

            gilCoords = (corot-np.min(lpoints), 2*np.pi/3, thetaGrid[-i])
            gxli, gyli, gzli = starToPlanet(gilCoords, planetCoords)
            smallestlFlow = flowLine3D(gxli, gyli, gzli, coordX, coordY, coordZ, error, vel, flowStep)

            s = s + (" " + str(np.min(lpoints)) + " " + str(np.min(smallestlFlow[2])) + " " + str(np.max(lpoints)) + " " + str(np.min(largestlFlow[2])))
        if(s == str(z)):
            f.write(str(z) + " 0 -1 0 -1")
        else:
           f.write(s)