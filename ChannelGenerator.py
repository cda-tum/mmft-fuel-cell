from HelperFunctions import *
import numpy as np
from scipy.spatial import Delaunay
import triangle as tr

def main():
    #Channels_2D()
    Channels_3D()
    #dim = [24e-7,24e-7,5e-7]
    #testchannel = incidentCubes(dim, 0.25*np.pi, 3e-7, 6e-7, 3e-7)
    #testchannel.fitGeometry()
    #testchannel.triangulate()
    #filename = "channels3D/Atestchannel.stl"
    #testchannel.write(filename)

def Channels_2D():
    origin = [0.0,0.0,0.0]
    dim = [5.0e-8,16.0e-6,2.0e-6]
    Heights = [25,50,75,100,125,150]
    Widths = [25,50,75,100]
    Frequencies = [1,2,4,8,16]
    Types = ["Wedge","Ramp","ReverseRamp","Step"]

    for type in Types:
        for frequency in Frequencies:
            for height in Heights:
                for width in Widths:
                    stepDim = [5.0e-8,width*1e-8,height*1e-8]
                    channel = Channel(type,dim,stepDim,frequency)
                    filename = "channels/{}_Freq{}_H{}_W{}.stl".format(type,frequency,height,width)
                    channel.generateChannel()
                    channel.write(filename)

def Channels_3D():
    origin = [0.0,0.0,0.0]
    channelWidths = [0]
    channelLengths = [30,60,120]
    Heights = [20,30,40]
    Lengths = [30,60,120]
    Widths = [30,60,120]
    for channelWidth in channelWidths:
        for channelLength in channelLengths:
            for height in Heights:
                for length in Lengths:
                    for width in Widths:
                        dim = [(channelWidth+width)*1e-8,(channelLength+length)*1e-8,0.5e-6]
                        stepDim = [width*1e-8,length*1e-8,height*1e-8]
                        cube = OffsetCube(origin,dim,stepDim)
                        #filename = "channels3D/{}_cW{}_cL{}_H{}_L{}_W{}.stl".format("Cube",\
                        #    channelWidth,channelLength,height,length,width)
                        #filename = "channels3D/{}_cL{}_H{}_L{}.stl".format("Step",\
                        #    channelLength,height,length)
                        #filename = "channels3D/{}_cW{}_H{}_W{}.stl".format("Groove",\
                        #    channelWidth,height,width)
                        filename = "channels3D/{}_cW{}_cL{}_H{}_L{}_W{}.stl".format("Offset",\
                            channelWidth,channelLength,height,length,width)

                        cube.write(filename)


if __name__ == "__main__":
    main()
