from HelperFunctions import *
import numpy as np
from scipy.spatial import Delaunay
import triangle as tr

def main():
    Channels_3D()

def Channels_3D():
    """3D channels with a selection of nanostructures: Cubes, Steps, Grooves, and Offset Cubes.
    The nanostructure can be set in Type."""

    # The four different types.
    Types = ["Cube","Step","Groove","Offset"]
    # Select here the nanostructure type to be generated.
    Type = Types[0]

    origin = [0.0,0.0,0.0]

    # Variables to for defining nanostructure dimension and frequency.
    channelWidths = [30,60,120]
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
                        match Type:
                            case "Cube":
                                cube = Cube(origin,dim,stepDim)
                                filename = "ChannelGeometries/{}_cW{}_cL{}_H{}_L{}_W{}.stl".format("Cube",\
                                channelWidth,channelLength,height,length,width)
                            case "Step":
                                cube = Cube(origin,dim,stepDim)
                                filename = "ChannelGeometries/{}_cL{}_H{}_L{}.stl".format("Step",\
                                channelLength,height,length)
                            case "Groove":
                                cube = Cube(origin,dim,stepDim)
                                filename = "ChannelGeometries/{}_cW{}_H{}_W{}.stl".format("Groove",\
                                channelWidth,height,width)
                            case "Offset":
                                cube = OffsetCube(origin,dim,stepDim)
                                filename = "ChannelGeometries/{}_cW{}_cL{}_H{}_L{}_W{}.stl".format("Offset",\
                                    channelWidth,channelLength,height,length,width)
                        cube.write(filename)


if __name__ == "__main__":
    main()
