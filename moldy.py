import pandas as pd #data management
import cv2 #image and video
import matplotlib.pyplot as plt #creating images (rendering)
import os #file access
import time #reporting code timings
from math import sin, cos
import numpy as np

VMD_PATH = "/Applications/VMD 2.0.0a7-pre2.app/Contents/vmd2/lib/vmd_MACOSXARM64"
MOLDY_PATH = ""

def renderBeads(beadOutput="beadOutput.txt",noFrames=600):    
    """
    
    A simple display function. Reads the bead output from moldy and renders
    images of beads and springs from the simulation using matplotlib

    Parameters
    ----------
    beadOutput : filepath, optional
        Should be left unchanged unless the user is providing their own bead 
        output file. The default is "beadOutput.txt".

    Returns
    -------
    imgArr : array
        Array of images for passing to the video creator function.

    """
    timeStart = time.time()
    print("Generating image frames...")
    imgArr = []
    outputDataBeads = pd.read_csv("beadOutput.txt", sep=' ')
    outputDataSprings = pd.read_csv("springOutput.txt", sep=' ')
    outputData1 = outputDataBeads.loc[outputDataBeads["ID"]==0]
    noBeads = outputDataBeads["ID"].max() + 1
    noSprings = outputDataSprings["ID"].max() + 1
    
    maxForce = outputDataSprings["Force"].max()
    # if maxForce == 0:
    #     maxForce = 1 #prevent divide by 0 errors
    
    Xbounds = [outputDataBeads["X"].min(), outputDataBeads["X"].max()]
    Ybounds = [outputDataBeads["Y"].min(), outputDataBeads["Y"].max()]
    Zbounds = [outputDataBeads["Z"].min(), outputDataBeads["Z"].max()]
    maxBeadRadius = outputDataBeads["Radius"].max()
    overallMin = min(Xbounds[0],Ybounds[0],Zbounds[0]) - maxBeadRadius
    overallMax = max(Xbounds[1],Ybounds[1],Zbounds[1]) + maxBeadRadius
    
    colours = ['brown','red',"gray",'b']
    
    f = 0
    
    
    times = outputData1["Time"].tolist()
    length = len(times)
    mult = int(length/noFrames)
    for t in times:
        if f%mult==0:
            #print(f)
            #print(t)
            fig, (XY, XZ) = plt.subplots(2,sharex=True,figsize=(5, 10))
            for b in range(0,noBeads):
                printData = outputDataBeads.loc[outputDataBeads["ID"]==b]
                if len(np.where(printData["Time"]==t)[0])!=0:
                    #print(b, t)
                    i = np.where(printData["Time"]==t)[0][0]
                else:
                    i = -1
                if i != -1:
                    X = printData["X"].iloc[i]
                    Y = printData["Y"].iloc[i]
                    Z = printData["Z"].iloc[i]
                    spec = printData["Species"].iloc[i]
                    #print(spec,t)
                    beadXY = plt.Circle([X,Y], 
                                      printData["Radius"].iloc[i], color=colours[spec],alpha=0.75,
                                      linewidth=1)
                    beadXY.set_ec("k")
                    beadXZ = plt.Circle([X,Z], 
                                      printData["Radius"].iloc[i], color=colours[spec],alpha=0.75,
                                      linewidth=1)
                    beadXZ.set_ec("k")
                    XY.add_patch(beadXY)
                    XZ.add_patch(beadXZ)
                
                
            XY.set_xlim([overallMin,overallMax])
            XY.set_ylim([overallMin,overallMax])
            XZ.set_xlim([overallMin,overallMax])
            XZ.set_ylim([overallMin,overallMax])
            
            XY.tick_params("both",direction='in')
            XY.tick_params("both",direction='in')
            
            XY.set_ylabel("Y (m)", fontsize=11)
            XZ.set_ylabel("Z (m)")
            XZ.set_xlabel("X (m)")
            XY.set_aspect('equal')
            XZ.set_aspect('equal')
            
            fig.tight_layout()
            #plt.show()
            plt.savefig("./graphicsOutput/images/frame"+str(i)+".png")
            imgArr.append(cv2.imread("./graphicsOutput/images/frame"+str(i)+".png"))
            plt.close()
            plt.clf()
        f += 1
    timeEnd = time.time()
    print("> %dm %.6fs" % ((timeEnd-timeStart)//60,(timeEnd-timeStart)%60))
    return imgArr

def renderBeadsSprings(beadOutput="beadOutput.txt",springOutput="springOutput.txt"):    
    """
    DEPRECATED DO NOT USE
    A simple display function. Reads the bead output from moldy and renders
    images of beads and springs from the simulation

    Parameters
    ----------
    beadOutput : filepath, optional
        Should be left unchanged unless the user is providing their own bead 
        output file. The default is "beadOutput.txt".
    springOutput : filepath, optional
        Should be left unchanged unless the user is providing their own bead 
        output file. The default is "springOutput.txt".

    Returns
    -------
    imgArr : array
        Array of images for passing to the video creator function.

    """
    timeStart = time.time()
    print("Generating image frames...")
    imgArr = []
    outputDataBeads = pd.read_csv("beadOutput.txt", sep=' ')
    outputDataSprings = pd.read_csv("springOutput.txt", sep=' ')
    outputData1 = outputDataBeads.loc[outputDataBeads["ID"]==0]
    noBeads = outputDataBeads["ID"].max() + 1
    noSprings = outputDataSprings["ID"].max() + 1
    
    maxForce = outputDataSprings["Force"].max()
    # if maxForce == 0:
    #     maxForce = 1 #prevent divide by 0 errors
    
    Xbounds = [outputDataBeads["X"].min(), outputDataBeads["X"].max()]
    Ybounds = [outputDataBeads["Y"].min(), outputDataBeads["Y"].max()]
    Zbounds = [outputDataBeads["Z"].min(), outputDataBeads["Z"].max()]
    maxBeadRadius = outputDataBeads["Radius"].max()
    overallMin = min(Xbounds[0],Ybounds[0],Zbounds[0]) - maxBeadRadius
    overallMax = max(Xbounds[1],Ybounds[1],Zbounds[1]) + maxBeadRadius
    
    colours = []
    for c in range(0,noBeads):
        colours.append([np.random.uniform(),np.random.uniform(),np.random.uniform()])
    
    for i in range(0,outputData1.shape[0]):
        fig, (XY, XZ) = plt.subplots(2,sharex=True,figsize=(5, 10))
        for b in range(0,noBeads):
            printData = outputDataBeads.loc[outputDataBeads["ID"]==b]
            X = printData["X"].iloc[i]
            Y = printData["Y"].iloc[i]
            Z = printData["Z"].iloc[i]
            beadXY = plt.Circle([X,Y], 
                              printData["Radius"].iloc[i], color=colours[b],alpha=0.75,
                              linewidth=1)
            beadXY.set_ec("k")
            beadXZ = plt.Circle([X,Z], 
                              printData["Radius"].iloc[i], color=colours[b],alpha=0.75,
                              linewidth=1)
            beadXZ.set_ec("k")
            XY.add_patch(beadXY)
            XZ.add_patch(beadXZ)
            
        for s in range(0,noSprings):
            printData = outputDataSprings.loc[outputDataSprings["ID"]==s]
            startInd = printData["Start"].iloc[i]
            endInd = printData["End"].iloc[i]
            
            beadStartData = outputDataBeads.loc[outputDataBeads["ID"]==startInd]
            beadEndData = outputDataBeads.loc[outputDataBeads["ID"]==endInd]
            
            startX = beadStartData["X"].iloc[i]
            startY = beadStartData["Y"].iloc[i]
            startZ = beadStartData["Z"].iloc[i]
            
            endX = beadEndData["X"].iloc[i]
            endY = beadEndData["Y"].iloc[i]
            endZ = beadEndData["Z"].iloc[i]
            
            forceOverMax = printData["Force"].iloc[i]/maxForce
            
            springXY = plt.Line2D([startX,endX],[startY,endY], 
                                  color=[forceOverMax,0,0],
                                  linewidth=1,
                                  marker='o')
            springXZ = plt.Line2D([startX,endX],[startZ,endZ], 
                                  color=[forceOverMax,0,0],
                                  linewidth=1,
                                  marker='o')
            springXY.set_zorder(1)
            springXZ.set_zorder(1)
            XY.add_line(springXY)
            XZ.add_line(springXZ)
            
            
        XY.set_xlim([overallMin,overallMax])
        XY.set_ylim([overallMin,overallMax])
        XZ.set_xlim([overallMin,overallMax])
        XZ.set_ylim([overallMin,overallMax])
        
        XY.tick_params("both",direction='in')
        XY.tick_params("both",direction='in')
        
        XY.set_ylabel("Y (m)", fontsize=11)
        XZ.set_ylabel("Z (m)")
        XZ.set_xlabel("X (m)")
        XY.set_aspect('equal')
        XZ.set_aspect('equal')
        
        fig.tight_layout()
        #plt.show()
        plt.savefig("./graphicsOutput/images/frame"+str(i)+".png")
        imgArr.append(cv2.imread("./graphicsOutput/images/frame"+str(i)+".png"))
        plt.close()
        plt.clf()
    timeEnd = time.time()
    print("> %dm %.6fs" % ((timeEnd-timeStart)//60,(timeEnd-timeStart)%60))
    return imgArr
   
def imagesToVideo(imageArray, videoName = "video.mp4"):
    """
    .mp4 video creator, taking image arrays as input. Video is created in the
    graphicsOutput folder.

    Parameters
    ----------
    imageArray : array
        An array of images as created by the render...() functions.
    videoName : string, optional
        final name of video file, extension must be mp4. 
        The default is "video.mp4".

    Returns
    -------
    None.

    """
    timeStart = time.time()
    print("Creating video from images...")
    height,width,layers=imageArray[0].shape
    
    fourcc = cv2.VideoWriter_fourcc(*'mp4v') 
    video = cv2.VideoWriter('./graphicsOutput/'+videoName,
                            fourcc, 30, (width, height))
    
    for img in imageArray:
        video.write(img)
    cv2.destroyAllWindows()
    video.release()
    timeEnd = time.time()
    print("output: " + "./graphicsOutput/"+videoName)
    print("> %dm %.6fs" % ((timeEnd-timeStart)//60,(timeEnd-timeStart)%60))
    
def runSimulation():
    """
    Runs the moldy executable from python.

    Returns
    -------
    None.

    """
    timeStart = time.time()
    print("Running moldy...")
    os.system("./moldy")
    timeEnd = time.time()
    print("> %dm %.6fs" % ((timeEnd-timeStart)//60,(timeEnd-timeStart)%60))

def createBead(units="SI",position=[0,0,0],velocity=[0,0,0],mass=0,radius=0,charge=0,species=0):
    """
    Create a single bead with no connections.

    Parameters
    ----------
    units : string, optional
        Input units (only SI is currently supported). The default is "SI".
    position : array, optional
        bead starting position. The default is [0,0,0].
    velocity : array, optional
        bead starting velocity. The default is [0,0,0].
    mass : float, optional
        bead mass. The default is 0.
    radius : float, optional
        bead radius. The default is 0.
    charge : float, optional
        bead charge. The default is 0.

    Returns
    -------
    str
        string formatted correctly for sending to the initialisation file.

    """
    if units == "SI":
        bead = ("%e %e %e %e %e %e %e %e %e %d" 
                         %(position[0],position[1],position[2],
                           velocity[0],velocity[1],velocity[2],
                           mass, radius, charge, species))
    return bead+'\n'

def createSpring(units="SI",startInd=0,endInd=1,natLength=0,springConst=0):
    """
    Create a spring between two beads by index.

    Parameters
    ----------
    units : string, optional
        Input units (only SI is currently supported). The default is "SI".
    startInd : int, optional
        index of the bead at one end of the spring. The default is 0.
    endInd : int, optional
        index of the bead at one end of the spring. The default is 0.
    natLength : float, optional
        Natural length of the spring. The default is 0.
    springConst : float, optional
        Hookean spring constant (spring stiffness). The default is 0.

    Returns
    -------
    str
        string formatted correctly for sending to the initialisation file.

    """
    if units == "SI":
        spring = ("%d %d %e %e" % (startInd,endInd,natLength,springConst))
    return spring+'\n'

def createHomoChain(n,units="SI",position=[0,0,0],velocity=[0,0,0],
                    mass=0,radius=0,charge=0,natLength=0,springConst=0,
                    chainStartInd=0,randWalk=False,species=0):
    """
    
    Creates a chain of homogenous beads and spring (all beads have identical 
    qualities, all springs have identical qualities).

    Parameters
    ----------
    n : int
        Number of monomers (chain length).
    units : string, optional
        Input units (only SI is currently supported). The default is "SI".
    position : array, optional
        bead starting position. The default is [0,0,0].
    velocity : array, optional
        bead starting velocity. The default is [0,0,0].
    mass : float, optional
        bead mass. The default is 0.
    radius : float, optional
        bead radius. The default is 0.
    charge : float, optional
        bead charge. The default is 0.
    natLength : float, optional
        Natural length of the spring. The default is 0.
    springConst : float, optional
        Hookean spring constant (spring stiffness). The default is 0.

    Returns
    -------
    beadArr : arr
        array of formatted strings for sending to the simulation initialisation
    springArr : arr
        array of formatted strings for sending to the simulation initialisation

    """
    beadArr = []
    springArr = []
    if units=="SI":
        for i in range(0,n):
            if randWalk:
                polar = np.random.uniform(-np.pi*0.5,+np.pi*0.5)
                azith = np.random.uniform(0,+np.pi*2)
                beadPosition = [natLength*sin(azith)*cos(polar),
                                natLength*sin(azith)*cos(azith),
                                natLength*cos(polar)]
                beadArr.append(createBead(position=[position[0]+beadPosition[0],
                                                    position[1]+beadPosition[1],
                                                    position[2]+beadPosition[2]], 
                                          velocity=[0,0,0],
                                          mass=mass,
                                          radius=radius,
                                          charge=charge,
                                          species=species))
                
            else:
                beadArr.append(createBead(position=[position[0]+i*natLength,
                                                    position[1],
                                                    position[2]],
                                          velocity=[0,0,0],
                                          mass=mass,
                                          radius=radius,
                                          charge=charge,
                                          species=species))
            if i != 0:
                springArr.append(createSpring(startInd=i-1+chainStartInd,
                                              endInd=i+chainStartInd,
                                              natLength=natLength,
                                              springConst=springConst))
        
    
    return beadArr, springArr

def createHomoLoop(n,units="SI",position=[0,0,0],velocity=[0,0,0],
                    mass=0,radius=0,charge=0,natLength=0,springConst=0):
    """
    
    Creates a loop of homogenous beads and spring (all beads have identical 
    qualities, all springs have identical qualities). I.E. polymer ends are 
    connected.

    Parameters
    ----------
    n : int
        Number of monomers (chain length).
    units : string, optional
        Input units (only SI is currently supported). The default is "SI".
    position : array, optional
        bead starting position. The default is [0,0,0].
    velocity : array, optional
        bead starting velocity. The default is [0,0,0].
    mass : float, optional
        bead mass. The default is 0.
    radius : float, optional
        bead radius. The default is 0.
    charge : float, optional
        bead charge. The default is 0.
    natLength : float, optional
        Natural length of the spring. The default is 0.
    springConst : float, optional
        Hookean spring constant (spring stiffness). The default is 0.

    Returns
    -------
    beadArr : arr
        array of formatted strings for sending to the simulation initialisation
    springArr : arr
        array of formatted strings for sending to the simulation initialisation

    """
    beadArr = []
    springArr = []
    d_theta = 2*3.14159/n
    dist = radius/sin(d_theta/2)
    if units=="SI":
        for i in range(0,n):
            beadArr.append(createBead(position=[dist*sin(i*d_theta),dist*cos(i*d_theta),0],
                                      velocity=[0,0,0],
                                      mass=mass,
                                      radius=radius,
                                      charge=charge))
            if i != 0:
                springArr.append(createSpring(startInd=i-1,
                                              endInd=i,
                                              natLength=natLength,
                                              springConst=springConst))
            else:
                springArr.append(createSpring(startInd=i,
                                              endInd=n-1,
                                              natLength=natLength,
                                              springConst=springConst))
        
    
    return beadArr, springArr
        
def initialiseObjects(beadArray, springArray, io='w'):
    """
    Sends beads and springs to the simulation files prior to running the
    simulation.

    Parameters
    ----------
    beadArray : array
        array of beads, generated from one of the creator files.
    springArray : array
        array of beads, generated from one of the creator files.
    io : char, optional
        method with which to open the initialisation files (clearing before
        writing or just appending). The default is 'w'.

    Returns
    -------
    None.

    """
    with open("./beadInit.txt", io) as f:
        for b in beadArray:
            f.write(b)
    with open("./springInit.txt", io) as f:
        for s in springArray:
            f.write(s)
            
def setSimulationParameters(maxTime,timestep,viscosity,temperature,
                            overlapStrength):
    """
    Set simulation parameters prior to running the simulation. Only required to
    be called when simulation parameters are changed.

    Parameters
    ----------
    maxTime : float
        maximum time in simulated seconds the simulation will run up to.
    timestep : float
        dt in the simulation.
    viscosity : float
        viscosity of the medium in which the beads and springs exist.
    temperature : float
        simulation temperature.
    overlapStrength : float
        strength with which beads resist overlapping. User defined, high values
        can affect the stability of the simulation.

    Returns
    -------
    None.

    """
    with open("./simulationParameters.txt", 'w') as f:
        f.write("%e %e %e %e %e" % (maxTime,timestep,viscosity,temperature,
                                       overlapStrength))
        
def peelLastFrame(writeBeadsFile, writeSpringFile, mass, springConst, natLength, charge=0, 
                  beadOutput="beadOutput.txt",springOutput="springOutput.txt"):
    outputDataBeads = pd.read_csv("beadOutput.txt", sep=' ')
    outputDataSprings = pd.read_csv("springOutput.txt", sep=' ')
    outputData1 = outputDataBeads.loc[outputDataBeads["ID"]==0]
    times = outputData1["Time"].tolist()
    
    beadsFinalTime = outputDataBeads.loc[outputDataBeads["Time"]==times[-1]]    
    
    beadX = beadsFinalTime["X"].tolist()
    beadX = np.subtract(beadX,np.mean(beadX))
    beadY = beadsFinalTime["Y"].tolist()
    beadY = np.subtract(beadY,np.mean(beadY))
    beadZ = beadsFinalTime["Z"].tolist()
    beadZ = np.subtract(beadZ,np.mean(beadZ))
    
    beadMass = mass
    beadRadius = beadsFinalTime["Radius"].tolist()
    beadSpecies = beadsFinalTime["Species"].tolist()
    
    springsFinalTime = outputDataSprings.loc[outputDataSprings["Time"]==times[-1]]    
    
    startInd = springsFinalTime["Start"].tolist()
    endInd = springsFinalTime["End"].tolist()
    
    with open(writeBeadsFile, 'w') as f:
        for b in range(0,len(beadX)):
            bead = ("%e %e %e %e %e %e %e %e %e %d\n" 
                             %(beadX[b], beadY[b], beadZ[b],
                               0, 0, 0,
                               beadMass, beadRadius[b], charge, beadSpecies[b]))
            f.write(bead)
            
    with open(writeSpringFile, 'w') as f:
        for s in range(0,len(startInd)):
            spring = ("%d %d %e %e\n" % (startInd[s],endInd[s],natLength,springConst))
            f.write(spring)
    

# SCRIPTING DURING SIMULATION
def createBlankScript(dt, tMax):
    t = 0
    scriptArray = []
    while t < tMax:
        scriptArray.append("%e %s %d %s %s" % (t, "none", 0, "none", "none")) #time action number info
        t += dt
    return scriptArray

def addEvent(script, ind, action, number, info1, info2):
    script[ind] = ("%e %s %d %s %s" % (time, action, number, info1, info2))
    return script

def setSimulationScript(script):
    with open("./testWrites/simulationScript.txt", 'w') as f:
        for line in script:
            f.write(line+'\n')
            
# INTERACTION MATRIX
def createBlankInteractionMatrix(species):
    size = len(species)    
    interMatArr = np.zeros((size,size))
    return interMatArr, species

def editInteractionMatrix(interactionMatrix, interactionBetween, newEnergy):
    interMatArr = interactionMatrix[0]
    species = interactionMatrix[1]    
    #find location of interactors
    row = None
    column = None
    for s in range(0,len(species)):
        if species[s]==interactionBetween[0]:
            row = s
        if species[s]==interactionBetween[1]:
            column = s
            
    if row==None or column==None:
        print("Interactors listed are not part of the interaction matrix")
        return interactionMatrix
    else:
        interMatArr[row][column] = newEnergy
        return interMatArr, species
    
def symmetriseInteractionMatrix(interactionMatrix):
    interMatArr = interactionMatrix[0]
    shape = np.shape(interMatArr)
    for i in range(0,shape[0]):
        for j in range(0,shape[1]):
            if interMatArr[i][j]!=0:
                interMatArr[j][i] = interMatArr[i][j]
    return interMatArr, interactionMatrix[1]
    
def setInteractionMatrix(interactionMatrix):
    interMatArr = interactionMatrix[0]
    species = interactionMatrix[1]
    with open("./testWrites/interactionMatrixARR.txt", 'w') as f:
        f.write(' '.join(species))
        f.write('\n')
        for row in interMatArr:
            f.write(' '.join(row.astype(str))+'\n')
            
            
# MISC.
def speciesIndMap(interactionMatrix,speciesName):
    ind = False
    species = interactionMatrix[1]
    for s in range(0,len(species)):
        if species[s]==speciesName:
            ind = s
    if ind:
        print(ind)
        return ind
    else:
        print("ind not found")
        
def directoryToVideo(directoryPath,videoName):
    import cv2
    #Read in images
    imgArr = []
    cwdList = os.listdir(directoryPath)
    for im in np.sort(cwdList):
        if im != ".DS_Store":
            path = os.path.join(directoryPath, im)
            imgArr.append(cv2.imread(path))
    #Stitch into video
    height,width,layers=imgArr[0].shape
    fourcc = cv2.VideoWriter_fourcc(*'mp4v') 
    video = cv2.VideoWriter(os.path.join(directoryPath, videoName),
                            fourcc, 30, (width, height))
    for img in imgArr:
        #print(img)
        video.write(img)
    cv2.destroyAllWindows()
    video.release()
 
    

def outputToVTF(species, beadOutput="beadOutput.txt", springOutput="springOutput.txt", noFrames=200,
                launchVMD=False):
    fileLines = []
    outputDataBeads = pd.read_csv(beadOutput, sep=' ')
    outputDataSprings = pd.read_csv(springOutput, sep=' ')
    outputData1 = outputDataBeads.loc[outputDataBeads["ID"]==0]
    noBeads = outputDataBeads["ID"].max() + 1
    noSprings = outputDataSprings["ID"].max() + 1
    # if maxForce == 0:
    #     maxForce = 1 #prevent divide by 0 errors
    
    Xbounds = [outputDataBeads["X"].min()*10**10, outputDataBeads["X"].max()*10**10]
    Ybounds = [outputDataBeads["Y"].min()*10**10, outputDataBeads["Y"].max()*10**10]
    Zbounds = [outputDataBeads["Z"].min()*10**10, outputDataBeads["Z"].max()*10**10]
    maxBeadRadius = outputDataBeads["Radius"].max()*10**10
    overallMin = min(Xbounds[0],Ybounds[0],Zbounds[0]) - maxBeadRadius
    overallMax = max(Xbounds[1],Ybounds[1],Zbounds[1]) + maxBeadRadius
    
    #bounding box
    boxLines = []
    boxLines.append("atom 0:7 radius 0.0 segid box\n")
    boxLines.append("bond 0:3,0::3\n")
    boxLines.append("bond 4:7,4::7\n")
    boxLines.append("bond 0:4,1:5,2:6,3:7\n")
    boxLines.append("timestep ordered\n")
    
    #create "atoms"
    fileLines.append("# STRUCTURE BLOCK\n")
    for i in range(0,noBeads):
        beadInfo = (outputDataBeads.loc[outputDataBeads["ID"]==i]).iloc[0]
        beadRadius = beadInfo["Radius"] * 10**10 #VTF radii are in angstroms
        beadSpecies = species[int(beadInfo["Species"])]
        currLine = "atom " + str(i)
        currLine += " radius" + (" %.2f" % beadRadius)
        currLine += " segid " + beadSpecies + '\n'
        fileLines.append(currLine)
    #create "bonds"    
    for j in range(0,noSprings):
        springInfo = (outputDataSprings.loc[outputDataBeads["ID"]==j]).iloc[0]
        springStart = int(springInfo["Start"])
        springEnd = int(springInfo["End"])
        currLine = "bond " + str(springStart) + ':' + str(springEnd) + '\n'
        fileLines.append(currLine)
    fileLines.append('\n')
        
    #populate with timestep data
    fileLines.append("# TIMESTEP BLOCKS\n")
    f=0
    times = outputData1["Time"].tolist()
    length = len(times)
    mult = int(length/noFrames)
    for t in times:
        if f%mult==0:
            #print(f)
            fileLines.append("# t = %e seconds\n" % t)
            fileLines.append("timestep indexed\n")
            idList = range(0,noBeads)
            #print(t)
            currData = outputDataBeads.loc[outputDataBeads["Time"]==t]
            #print(currData)
            #print(idList)
            removedInd = 0
            for b in range(0,noBeads):
                currBead = currData.loc[currData["ID"]==b]
                if not currBead.empty:
                    X = currBead["X"].iloc[0] * 10**10
                    Y = currBead["Y"].iloc[0] * 10**10
                    Z = currBead["Z"].iloc[0] * 10**10
                    fileLines.append(str(b) + ' '
                                    +str(X) + ' '
                                    +str(Y) + ' '
                                    +str(Z) + ' ' + '\n')
                    idList = np.delete(idList,removedInd)
                    #print(idList)
                else:
                    #preparation for dealing with unspawned beads...
                    fileLines.append(str(b) + ' '
                                    +str(-1e9) + ' '
                                    +str(-1e9) + ' '
                                    +str(-1e9) + ' ' + '\n')
                #print(currBead)
            
            boxLines.append("timestep ordered\n")
            boxLines.append(str(Xbounds[0])+' '+str(Ybounds[0])+' '+str(Zbounds[0])+'\n')
            boxLines.append(str(Xbounds[1])+' '+str(Ybounds[0])+' '+str(Zbounds[0])+'\n')
            boxLines.append(str(Xbounds[1])+' '+str(Ybounds[1])+' '+str(Zbounds[0])+'\n')
            boxLines.append(str(Xbounds[0])+' '+str(Ybounds[1])+' '+str(Zbounds[0])+'\n')
            boxLines.append(str(Xbounds[0])+' '+str(Ybounds[0])+' '+str(Zbounds[1])+'\n')
            boxLines.append(str(Xbounds[1])+' '+str(Ybounds[0])+' '+str(Zbounds[1])+'\n')
            boxLines.append(str(Xbounds[1])+' '+str(Ybounds[1])+' '+str(Zbounds[1])+'\n')
            boxLines.append(str(Xbounds[0])+' '+str(Ybounds[1])+' '+str(Zbounds[1])+'\n')
        f += 1
                
    with open("./testWrites/molecule_output.vtf", 'w') as f:
        for l in fileLines:
            f.write(l)
            
    
    with open("./testWrites/boundingBox.vtf", 'w') as f:
        for l in boxLines:
            f.write(l)
    
    boundingBoxLoc = os.path.realpath("./testWrites/boundingBox.vtf")
    moleculeLoc = os.path.realpath("./testWrites/molecule_output.vtf")
    #generate tcl script
    tclLines = []
    tclLines.append("display resetview\n")
    tclLines.append("color Axes Labels gray\n")
    tclLines.append("display projection Perspective\n")
    
    
    tclLines.append("color Display Background white\n")
    tclLines.append("display backgroundgradient off\n")
    tclLines.append("display depthcue off\n")
    tclLines.append("mol new %s type {vtf} first 0 last -1 step 1 waitfor -1\n" %
                    boundingBoxLoc)
    tclLines.append("mol new %s type {vtf} first 0 last -1 step 1 waitfor -1\n" %
                    moleculeLoc)
    tclLines.append("color Segname box black\n")
    tclLines.append("mol modcolor 0 0 SegName\n")
    tclLines.append("mol modcolor 0 1 SegName\n")
    
    tclLines.append("mol modstyle 0 1 VDW 1.000000 16.000000\n")
    tclLines.append("mol modmaterial 0 1 Diffuse\n")
    
    tclLines.append("rotate x by -90\n")
    tclLines.append("rotate y by -135\n")    
    tclLines.append("rotate x by 45\n")
    tclLines.append("wait 18000") #required to stop it closing instantly
    
    with open("./testWrites/moldyOutputSetup.tcl", 'w') as f:
        for l in tclLines:
            f.write(l)
            
    if launchVMD:
        assert(len(VMD_PATH)!=0), "VMD_PATH has length zero: Make sure VMD_PATH points to the VMD executable."
        import subprocess
        molpath = os.path.realpath("./testWrites/moldyOutputSetup.tcl")
        command1 = f'"{VMD_PATH}" -e {molpath}'
        commands = [VMD_PATH,"-e",molpath]
        print(command1)
        p1 = subprocess.Popen(commands)
        print(p1)
    else:
        print("\n!!! IMPORTANT !!!")
        print("run the following command in VMD console:")
        print("source %s" % os.path.realpath("./testWrites/moldyOutputSetup.tcl"))            
    
# # EXAMPLE USAGE:
# # Below is an example set up. A loop is created of positive beads and a negative
# # bead is placed in the centre.
# # Simulation parameters are also set (in S.I.)
# # The simulation is run and the results are rendered into images. The images are
# # then stitched into a video.

# beads, springs = createHomoLoop(32,mass=1e-9,radius=1e-6,charge=1e-6,
#                                  natLength=1.9e-6, springConst=1e-7)
# negBead = createBead(position=[0,0,0],mass=1e-8,radius=1e-6,charge=-1e-2)

# beads = beads + [negBead]

setSimulationParameters(1e-6,1e-10,1e-3,300,1e-10)

# initialiseObjects(beads,springs)

# nChains = 3
# chainLength = 2
# beads = []
# springs = []
# spread = 1e-8
# for i in range(0,nChains):
#     beadPosition = [np.random.uniform(-spread,+spread),
#                     np.random.uniform(-spread,+spread),
#                     np.random.uniform(-spread,+spread)]
#     bead, spring = createHomoChain(chainLength,position=beadPosition,mass=1e-12,radius=1e-9,
#                            natLength=(2.0**(1/6))*2e-9,springConst=1e-5,chainStartInd=chainLength*i,
#                            randWalk=False)
#     beads = beads + bead
#     springs = springs + spring

# initialiseObjects(beads,springs)

#Sle 
beads = []
springs = []
springLength = (2.0**(1/6))*2e-9
#springLength = 2e-9
for i in range(0,7):
        
    s = 2
    
    if i<3:
        s = 0 #weak
    if i==6:
        s = 1 #strong
    else:
        pass

    beads.append(createBead(position=[i*springLength,0,0],
                            mass = 1e-18, radius=1e-9,species=s))
    if i>0:
        springs.append(createSpring(startInd=i-1,endInd=i,
                                    natLength=springLength,
                                    springConst=1e-1))

    

# beads2, springs2 = createHomoChain(7,position=[0,6e-9,0],mass=1e-18,radius=1e-9,natLength=springLength,
#                                  springConst=1e-3,species=3,chainStartInd=7)
# beads2, springs2 = createHomoChain(7,position=[0,6e-9,0],mass=1e-18,radius=1e-9,natLength=springLength,
#                                  springConst=1e-1   ,species=3,chainStartInd=0,randWalk=False)

beadSolo = createBead(position=[0,6e-9,6e-9],mass=1e-18,radius=1e-9,species=3)

#initialiseObjects(beads+beads2,springs+springs2)
initialiseObjects(beads,springs)
#initialiseObjects(beads2,springs2)
runSimulation()
imageArray = renderBeads(noFrames=100)
imagesToVideo(imageArray)   




