"""These are helper functions used to define channel objects that are written
in STL format"""

import numpy as np
from decimal import Decimal
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import triangle as tr

class Channel():
    """A channel object with the following arguments:
    type - The type of element (smooth, step, cavity, ramp, wedge or slit)
    dim - The dimension of the channel in np.array([x,y,z])
    stepDim=np.zeros(3) - The dimension of the step in np.array([x,y,z])
    interval=0.0 - The distance along the y-axis between the elements
    gamma=0.0 - The incidence angle for elements that are turned out of x-plane
    theta=0.0 - The inclination angle of a ramp element"""

    def __init__(self, type, dim, stepDim, frequency):

        self.type = type
        self.dim = dim
        self.stepDim = stepDim
        self.frequency = frequency
        self.triangles = []
        self.rectangles = []
        self.elements = []
        self.subChannels = []

        typeList = ["Smooth","Step","Wedge","Ramp","ReverseRamp"]
        assert (type in typeList), "Invalid wall type for bottom wall."

        assert (dim[0] == stepDim[0]), "The widths of the channel and step are not equal."

        if (type != "Smooth"):
            assert ((stepDim>np.zeros(3)).all), "Invalid step dimension."
            assert (frequency>=1 and frequency<=16), "Invalid element frequency."

    def generateChannel(self):
        # Front and Back walls
        origin = np.array( [[0,0,0], \
                            [0,self.dim[1],0]])
        size = np.array([[self.dim[0],self.dim[2]], \
                        [self.dim[2],self.dim[0]]])
        t1 = np.array([[1,0,0], \
                       [0,0,1]])
        t2 = np.array([[0,0,1], \
                       [1,0,0]])

        for x in range(2):
            self.rectangles.append(Rectangle(origin[x], size[x], t1[x], t2[x]))

        if (type == "Smooth"):
            self.subChannels.append(SubChannel(origin[0],self.dim))
        else:
            subChannelLength = (self.dim[1]/self.frequency)-self.stepDim[1]
            if (self.dim[1]-self.frequency*self.stepDim[1]) < -1e-12:
                subChannelLength = -1.0
            assert (subChannelLength >= 0.0), "Error: The stepDim*frequency is bigger than total channel length."
            # Add initial 1/2 subChannel
            self.subChannels.append(SubChannel(origin[0],np.array([self.dim[0],subChannelLength/2,self.dim[2]]),self.stepDim))
            # Add n-1 elements and subChannels
            for x in range(self.frequency-1):
                loc = subChannelLength/2 + x*(self.stepDim[1]+subChannelLength)
                match self.type:
                    case "Step":
                        self.elements.append(StepElement(np.array([0,loc,0]),self.dim,self.stepDim))
                    case "Wedge":
                        self.elements.append(WedgeElement(np.array([0,loc,0]),self.dim,self.stepDim))
                    case "Ramp":
                        self.elements.append(RampElement(np.array([0,loc,0]),self.dim,self.stepDim))
                    case "ReverseRamp":
                        self.elements.append(ReverseRampElement(np.array([0,loc,0]),self.dim,self.stepDim))
                self.subChannels.append(SubChannel(np.array([0,loc+self.stepDim[1],0]),np.array([self.dim[0],subChannelLength,self.dim[2]]),self.stepDim))
            # Add n-th element
            loc = subChannelLength/2 + (self.frequency-1)*(self.stepDim[1]+subChannelLength)
            match self.type:
                case "Step":
                    self.elements.append(StepElement(np.array([0,loc,0]),self.dim,self.stepDim))
                case "Wedge":
                    self.elements.append(WedgeElement(np.array([0,loc,0]),self.dim,self.stepDim))
                case "Ramp":
                    self.elements.append(RampElement(np.array([0,loc,0]),self.dim,self.stepDim))
                case "ReverseRamp":
                    self.elements.append(ReverseRampElement(np.array([0,loc,0]),self.dim,self.stepDim))
            # Add last 1/2 subChannel
            self.subChannels.append(SubChannel(np.array([0,loc+self.stepDim[1],0]),np.array([self.dim[0],subChannelLength/2,self.dim[2]]),self.stepDim))

    def write(self,filename):
        file = open(filename, "w")
        file.write("solid\n")
        for subChannel in self.subChannels:
            subChannel.write(file)
        for element in self.elements:
            element.write(file)
        for rectangle in self.rectangles:
            rectangle.write(file)
        for triangle in self.triangles:
            triangle.write(file)
        file.write("endsolid")
        file.close()

class Cube():
    """A cube element object, used for 3D channels with cubic obstacles"""
    def __init__(self, origin, Dim, stepDim):
        self.cubeOrigin = np.array([(Dim[0]-stepDim[0])/2, \
                                    (Dim[1]-stepDim[1])/2, \
                                    0.0])
        self.origin = origin
        self.dim = Dim
        self.stepDim = stepDim
        self.triangles = []
        self.rectangles = []

        # Vertex order, first big cube, then small cube
        #     5------6      13-----14
        #    /|     /|     /|     /|
        #   4------7 |    12-----15|
        #   | 1----|-2    | 9----|-10
        #   |/     |/     |/     |/
        #   0------3      8------11

        self.vertices = np.array([\
                                [0,0,0],\
                                [0,Dim[1],0],\
                                [Dim[0],Dim[1],0],\
                                [Dim[0],0,0],\
                                [0,0,Dim[2]],\
                                [0,Dim[1],Dim[2]],\
                                [Dim[0],Dim[1],Dim[2]],\
                                [Dim[0],0,Dim[2]],\
                                [self.cubeOrigin[0],self.cubeOrigin[1],0],\
                                [self.cubeOrigin[0],self.cubeOrigin[1]+stepDim[1],0],\
                                [self.cubeOrigin[0]+stepDim[0],self.cubeOrigin[1]+stepDim[1],0],\
                                [self.cubeOrigin[0]+stepDim[0],self.cubeOrigin[1],0],\
                                [self.cubeOrigin[0],self.cubeOrigin[1],stepDim[2]],\
                                [self.cubeOrigin[0],self.cubeOrigin[1]+stepDim[1],stepDim[2]],\
                                [self.cubeOrigin[0]+stepDim[0],self.cubeOrigin[1]+stepDim[1],stepDim[2]],\
                                [self.cubeOrigin[0]+stepDim[0],self.cubeOrigin[1],stepDim[2]]])

        for vertex in self.vertices:
            vertex += self.origin

        # Channel Front
        self.rectangles.append(Rectangle(self.vertices[0], np.array([Dim[0],Dim[2]]), np.array([1,0,0]), np.array([0,0,1])))
        # Channel Back
        self.rectangles.append(Rectangle(self.vertices[1], np.array([Dim[2],Dim[0]]), np.array([0,0,1]), np.array([1,0,0])))
        # Channel Left
        self.rectangles.append(Rectangle(self.vertices[0], np.array([Dim[2],Dim[1]]), np.array([0,0,1]), np.array([0,1,0])))
        # Channel Right
        self.rectangles.append(Rectangle(self.vertices[3], np.array([Dim[1],Dim[2]]), np.array([0,1,0]), np.array([0,0,1])))
        # Channel Top
        self.rectangles.append(Rectangle(self.vertices[4], np.array([Dim[0],Dim[1]]), np.array([1,0,0]), np.array([0,1,0])))
        # Cube Front
        self.rectangles.append(Rectangle(self.vertices[8], np.array([stepDim[2],stepDim[0]]), np.array([0,0,1]), np.array([1,0,0])))
        # Cube Back
        self.rectangles.append(Rectangle(self.vertices[13], np.array([stepDim[2],stepDim[0]]), np.array([0,0,-1]), np.array([1,0,0])))
        # Cube Left
        self.rectangles.append(Rectangle(self.vertices[8], np.array([stepDim[1],stepDim[2]]), np.array([0,1,0]), np.array([0,0,1])))
        # Cube Right
        self.rectangles.append(Rectangle(self.vertices[11], np.array([stepDim[2],stepDim[1]]), np.array([0,0,1]), np.array([0,1,0])))
        # Cube Top
        self.rectangles.append(Rectangle(self.vertices[12], np.array([stepDim[1],stepDim[0]]), np.array([0,1,0]), np.array([1,0,0])))
        # Channel Bottom
        self.triangles.append(Triangle(np.array([self.vertices[0],self.vertices[11],self.vertices[3]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[0],self.vertices[8],self.vertices[11]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[1],self.vertices[8],self.vertices[0]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[1],self.vertices[9],self.vertices[8]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[2],self.vertices[9],self.vertices[1]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[2],self.vertices[10],self.vertices[9]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[3],self.vertices[10],self.vertices[2]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[3],self.vertices[11],self.vertices[10]]),np.array([0,0,-1])))

    def write(self,filename):
        file = open(filename, "w")
        file.write("solid\n")
        for rectangle in self.rectangles:
            rectangle.write(file)
        for triangle in self.triangles:
            triangle.write(file)
        file.write("endsolid")
        file.close()

class OffsetCube():
        """An element of two cubes that are positions with an offset to each other"""
        def __init__(self, origin, Dim, cubeDim):
            self.origin = origin
            self.dim = Dim

            self.cubeElements = []
            self.rectangles = []
            self.triangles = []

            self.totalDim = np.array([2*Dim[0], 2*Dim[1], Dim[2]])

            self.cubeOrigin = np.array([(Dim[0]-cubeDim[0])/2, \
                                        (Dim[1]-cubeDim[1])/2, \
                                        0.0])

            offsetOrigin = self.origin + np.array([Dim[0], Dim[1], 0.0])

            self.vertices = np.array([\
                                    [0,0,0],\
                                    [0,Dim[1],0],\
                                    [0,2*Dim[1],0],\
                                    [Dim[0],0,0],\
                                    [Dim[0],Dim[1],0],\
                                    [Dim[0],2*Dim[1],0],\
                                    [2*Dim[0],0,0],\
                                    [2*Dim[0],Dim[1],0],\
                                    [2*Dim[0],2*Dim[1],0],\
                                    [0,0,Dim[2]],\
                                    [0,2*Dim[1],Dim[2]],\
                                    [2*Dim[0],2*Dim[1],Dim[2]],\
                                    [2*Dim[0],0,Dim[2]]])

            # Cubes
            self.cubeElements.append(CubeElement(self.origin, Dim, self.cubeOrigin, cubeDim))
            self.cubeElements.append(CubeElement(offsetOrigin, Dim, self.cubeOrigin, cubeDim))
            # Channel Bottom
            self.rectangles.append(Rectangle(self.vertices[3], np.array([Dim[1],Dim[0]]), np.array([0,1,0]), np.array([1,0,0])))
            self.rectangles.append(Rectangle(self.vertices[1], np.array([Dim[1],Dim[0]]), np.array([0,1,0]), np.array([1,0,0])))
            # Channel Front
            self.triangles.append(Triangle(np.array([self.vertices[6],self.vertices[12],self.vertices[3]]),np.array([0,-1,0])))
            self.triangles.append(Triangle(np.array([self.vertices[12],self.vertices[9],self.vertices[3]]),np.array([0,-1,0])))
            self.triangles.append(Triangle(np.array([self.vertices[3],self.vertices[9],self.vertices[0]]),np.array([0,-1,0])))
            # Channel Back
            self.triangles.append(Triangle(np.array([self.vertices[8],self.vertices[5],self.vertices[11]]),np.array([0,1,0])))
            self.triangles.append(Triangle(np.array([self.vertices[11],self.vertices[5],self.vertices[10]]),np.array([0,1,0])))
            self.triangles.append(Triangle(np.array([self.vertices[5],self.vertices[2],self.vertices[10]]),np.array([0,1,0])))
            # Channel Left
            self.triangles.append(Triangle(np.array([self.vertices[0],self.vertices[9],self.vertices[1]]),np.array([-1,0,0])))
            self.triangles.append(Triangle(np.array([self.vertices[9],self.vertices[10],self.vertices[1]]),np.array([-1,0,0])))
            self.triangles.append(Triangle(np.array([self.vertices[1],self.vertices[10],self.vertices[2]]),np.array([-1,0,0])))
            # Channel Right
            self.triangles.append(Triangle(np.array([self.vertices[6],self.vertices[7],self.vertices[12]]),np.array([1,0,0])))
            self.triangles.append(Triangle(np.array([self.vertices[7],self.vertices[11],self.vertices[12]]),np.array([1,0,0])))
            self.triangles.append(Triangle(np.array([self.vertices[7],self.vertices[8],self.vertices[11]]),np.array([1,0,0])))
            # Channel Top
            self.rectangles.append(Rectangle(self.vertices[9], np.array([2*Dim[0],2*Dim[1]]), np.array([1,0,0]), np.array([0,1,0])))

        def write(self,filename):
            file = open(filename, "w")
            file.write("solid\n")
            for cubeElement in self.cubeElements:
                cubeElement.write(file)
            for rectangle in self.rectangles:
                rectangle.write(file)
            for triangle in self.triangles:
                triangle.write(file)
            file.write("endsolid")
            file.close()

class Herringbone():
    """An element of a repetitive Herringbone structure"""
    #def __init__(self, origin, Dim, cubeDim):

    #def write(self, filename):

class incidentCubes():
    """A cube element with an angle of attack w.r.t. the incoming strea, used for Herringbone structure."""
    def __init__(self, dim, alpha, width, length, height):
        self.dim = dim
        self.alpha = alpha
        self.width = width
        self.length = length
        self.height = height

        self.vertices = {
            "Domain": [],
            "A": [],
            "B": [],
            "Front": [],
            "Left": [],
            "Back": [],
            "Right": [],
            "Top": []
        }

        self.triangles = []
        self.rectangles = []
        self.floorTriangles = []
        self.blockTriangles = []

        self.holes = np.array([])
        self.segments = np.array([])

        self.cut = np.array([0,0])

        self.domain = np.array([[0,0],
                                [0, dim[1]],
                                [dim[0], dim[1]],
                                [dim[0], 0]])
        self.segdomain = np.array([ [0,1],
                                    [1,2],
                                    [2,3],
                                    [3,0],])

        # Block a
        self.a = np.array([ [np.sin(alpha)*width, 0.0],
                            [0.0, np.cos(alpha)*width],
                            [np.cos(alpha)*length, np.sin(alpha)*length + np.cos(alpha)*width],
                            [np.cos(alpha)*length + np.sin(alpha)*width, np.sin(alpha)*length]
                            ])
        self.sega = self.segdomain

        # Block b
        self.b = np.array([ [np.cos(alpha)*length, 0.0],
                            [0.0, np.sin(alpha)*length],
                            [np.sin(alpha)*width, np.sin(alpha)*length + np.cos(alpha)*width],
                            [np.cos(alpha)*length + np.sin(alpha)*width, np.cos(alpha)*width]
                            ])
        self.segb = self.segdomain

        for vertex in self.b:
            vertex += np.array([0.5*dim[0], 0.5*dim[1]])

        xoffset = (dim[0] - self.b[3][0])/2
        yoffset = (dim[1] - self.b[2][1])/2

        for point in self.a:
            point += np.array([xoffset, yoffset])
        for point in self.b:
            point += np.array([xoffset, yoffset])

        self.holes = np.array([ [(self.a[1][0]+self.a[3][0])/2, (self.a[0][1]+self.a[2][1])/2],
                                [(self.b[1][0]+self.b[3][0])/2, (self.b[0][1]+self.b[2][1])/2]])

    def triangulate(self):

        pts = np.vstack([self.domain, self.a, self.b])
        seg = np.vstack([self.segdomain, self.sega + self.segdomain.shape[0],  self.segb + self.segdomain.shape[0] + self.sega.shape[0]])

        A = dict(vertices=pts, segments=seg, holes=self.holes)
        B = tr.triangulate(A, 'pA')

        nBottom = np.array([0,0,-1])
        nFront = np.array([0,-1,0])
        nLeft = np.array([-1,0,0])
        nBack = np.array([0,1,0])
        nRight = np.array([0,1,0])
        nTop = np.array([0,0,1])

        for triangle in B['triangles']:
            p1 = np.array([B['vertices'][triangle[0]][0], B['vertices'][triangle[0]][1], 0])
            p2 = np.array([B['vertices'][triangle[1]][0], B['vertices'][triangle[1]][1], 0])
            p3 = np.array([B['vertices'][triangle[2]][0], B['vertices'][triangle[2]][1], 0])
            self.floorTriangles.append(Triangle(np.array([p1,p3,p2]),nBottom))

        self.extrude()

        t = np.array([np.cos(self.alpha), np.sin(self.alpha), 0])
        n = np.array([np.sin(self.alpha), -np.cos(self.alpha), 0])

        # Block A is not cut
        if self.cut[0] == 0:
            self.addTriangle('A', [0,1,5], n)
            self.addTriangle('A', [5,4,0], n)
            self.addTriangle('A', [1,2,5], t)
            self.addTriangle('A', [2,6,5], t)
            self.addTriangle('A', [3,6,2], -n)
            self.addTriangle('A', [3,7,6], -n)
            self.addTriangle('A', [0,4,7], -t)
            self.addTriangle('A', [0,7,3], -t)
            self.addTriangle('A', [5,6,7], nBottom)
            self.addTriangle('A', [4,5,7], nBottom)

            self.vertices['Front'].append(self.vertices['Domain'][0])
            self.vertices['Front'].append(self.vertices['Domain'][4])
            self.vertices['Front'].append(self.vertices['Domain'][7])
            self.vertices['Front'].append(self.vertices['Domain'][3])
            #self.vertices['Front'].append(self.vertices['A'][0])
            #self.vertices['Front'].append(self.vertices['A'][4])

            self.addTriangle('Front', [0,1,2], nFront)
            self.addTriangle('Front', [2,3,0], nFront)

            self.vertices['Left'].append(self.vertices['Domain'][1])
            self.vertices['Left'].append(self.vertices['Domain'][5])
            self.vertices['Left'].append(self.vertices['Domain'][4])
            self.vertices['Left'].append(self.vertices['Domain'][0])
            #self.vertices['Left'].append(self.vertices['A'][1])
            #self.vertices['Left'].append(self.vertices['A'][5])

            self.addTriangle('Left', [0,1,2], nLeft)
            self.addTriangle('Left', [2,3,0], nLeft)

            self.vertices['Back'].append(self.vertices['Domain'][2])
            self.vertices['Back'].append(self.vertices['Domain'][6])
            self.vertices['Back'].append(self.vertices['Domain'][5])
            self.vertices['Back'].append(self.vertices['Domain'][1])

            self.addTriangle('Back', [0,1,2], nBack)
            self.addTriangle('Back', [2,3,0], nBack)

            self.vertices['Right'].append(self.vertices['Domain'][3])
            self.vertices['Right'].append(self.vertices['Domain'][7])
            self.vertices['Right'].append(self.vertices['Domain'][6])
            self.vertices['Right'].append(self.vertices['Domain'][2])

            self.addTriangle('Right', [0,1,2], nRight)
            self.addTriangle('Right', [2,3,0], nRight)

            self.vertices['Top'].append(self.vertices['Domain'][4])
            self.vertices['Top'].append(self.vertices['Domain'][5])
            self.vertices['Top'].append(self.vertices['Domain'][6])
            self.vertices['Top'].append(self.vertices['Domain'][7])

            self.addTriangle('Top', [0,1,2], nTop)
            self.addTriangle('Top', [2,3,0], nTop)

        # Block B is not cut
        if self.cut[1] == 0:
            self.addTriangle('B', [0,1,5], n)
            self.addTriangle('B', [5,4,0], n)
            self.addTriangle('B', [1,2,5], t)
            self.addTriangle('B', [2,6,5], t)
            self.addTriangle('B', [3,6,2], -n)
            self.addTriangle('B', [3,7,6], -n)
            self.addTriangle('B', [0,4,7], -t)
            self.addTriangle('B', [0,7,3], -t)
            self.addTriangle('B', [5,6,7], nBottom)
            self.addTriangle('B', [4,5,7], nBottom)

        # Block A is cut once
        #if self.cut[0] == 1:

        # Block B is cut once
        #if self.cut[0] == 1:

        # Block A is cut twice
        #if self.cut[0] == 1:

        # Block B is cut twice
        #if self.cut[0] == 1:

        tr.compare(plt, A, B)

        #plt.triplot(vertices[:,0], vertices[:,1], tri.simplices)
        #plt.plot(vertices[:,0], vertices[:,1], 'o')
        plt.show()

    def addTriangle(self, key, indices, vec):
        p1 = self.vertices[key][indices[0]]
        p2 = self.vertices[key][indices[1]]
        p3 = self.vertices[key][indices[2]]
        self.blockTriangles.append(Triangle(np.array([p1, p2, p3]), vec))

    def extrude(self):
        for vertex in self.domain:
            self.vertices["Domain"].append([vertex[0], vertex[1], 0])
        for vertex in self.domain:
            self.vertices["Domain"].append([vertex[0], vertex[1], self.dim[2]])
        for vertex in self.a:
            self.vertices["A"].append([vertex[0], vertex[1], 0])
        for vertex in self.a:
            self.vertices["A"].append([vertex[0], vertex[1], self.height])
        for vertex in self.b:
            self.vertices["B"].append([vertex[0], vertex[1], 0])
        for vertex in self.b:
            self.vertices["B"].append([vertex[0], vertex[1], self.height])

    def write(self,filename):
        file = open(filename, "w")
        file.write("solid\n")
        for rectangle in self.rectangles:
            rectangle.write(file)
        for triangle in self.floorTriangles:
            triangle.write(file)
        for triangle in self.blockTriangles:
            triangle.write(file)
        file.write("endsolid")
        file.close()

    def fitGeometry(self):
        for vertex in self.a:
            if vertex[1] > self.dim[1]:
                newValue = vertex[1]%self.dim[1]
                vertex[1] = newValue
                self.cut[0] += 1
        for vertex in self.b:
            if vertex[1] > self.dim[1]:
                newValue = vertex[1]%self.dim[1]
                vertex[1] = newValue
                self.cut[1] += 1

        # Block A is cut once
        if self.cut[0] == 1:
            excess = self.a[2][1]
            x1 = self.a[2][0] - excess/(np.tan(self.alpha))
            x2 = self.a[2][0] + excess/(np.tan(self.alpha))
            out = np.vstack((self.a, np.array([
                                        [x1, self.dim[1]],
                                        [x2, self.dim[1]],
                                        [x1, 0.0],
                                        [x2, 0.0]])))
            self.a = out
            self.sega = np.array([[0,1], [1,4], [4,5], [5,3], [3,0], [6,2], [2,7], [7,6]])
            holes = np.vstack([self.holes, np.array([(x1+x2)/2, excess/2])])
            self.holes = holes

        # Block A is cut twice
        if self.cut[0] == 2:
            excess = self.a[2][1]
            excess2 = self.a[3][1]
            x1 = self.a[2][0] - excess/(np.tan(self.alpha))
            x2 = self.a[2][0] + excess/(np.tan(self.alpha)) - 2*excess2/(np.tan(self.alpha))
            out = np.vstack((self.a, np.array([
                                        [x1, self.dim[1]],
                                        [x2, self.dim[1]],
                                        [x1, 0.0],
                                        [x2, 0.0]])))
            self.a = out
            self.sega = np.array([[0,1], [1,4], [4,5], [5,0], [6,2], [2,3], [3,7],
             [7,6]])
            holes = np.vstack([self.holes, np.array([[(x1+x2)/2,(self.width)/2]])])
            self.holes = holes

        # Block B is cut once
        if self.cut[1] == 1:
            excess = self.b[2][1]
            x1 = self.b[2][0] - excess/(np.tan(self.alpha))
            x2 = self.b[2][0] + excess/(np.tan(self.alpha))
            out = np.vstack((self.b, np.array([
                                        [x1, self.dim[1]],
                                        [x2, self.dim[1]],
                                        [x1, 0.0],
                                        [x2, 0.0]])))
            self.b = out
            self.segb = np.array([[0,1], [1,4], [4,5], [5,3], [6,2], [2,7], [7,6]])
            holes = np.vstack([self.holes, np.array([(x1+x2)/2, excess/2])])
            self.holes = holes

        # Block B is cut twice
        if self.cut[1] == 2:
            excess = self.b[2][1]
            excess2 = self.b[1][1]
            x1 = self.b[2][0] - excess/(np.tan(self.alpha)) + 2*(excess2/(np.tan(self.alpha)))
            x2 = self.b[2][0] + excess/(np.tan(self.alpha))
            out = np.vstack((self.b, np.array([
                                        [x1, self.dim[1]],
                                        [x2, self.dim[1]],
                                        [x1, 0.0],
                                        [x2, 0.0]])))
            self.b = out
            self.segb = np.array([[0,4], [4,5], [5,3], [3,0], [6,1], [1,2], [2,7], [7,6]])
            holes = np.vstack([self.holes, np.array([(x1+x2)/2,(self.width)/2])])
            self.holes = holes


class StepElement():
    """An step elememt object, used to describe the bottom wall pattern for channel step"""
    def __init__(self, origin, Dim, stepDim):
        self.rectangles = []
        self.triangles = []
        self.origin = np.array( [[0,0,stepDim[2]],\
                                [0,stepDim[1],stepDim[2]],\
                                [0,0,stepDim[2]],\
                                [stepDim[0],0,Dim[2]],\
                                [0,0,Dim[2]]])
        # Step front
        self.rectangles.append(Rectangle(origin, np.array([stepDim[2],stepDim[0]]), \
        np.array([0,0,1]), np.array([1,0,0])))
        # Step top
        self.rectangles.append(Rectangle(origin+self.origin[0], np.array([stepDim[1],stepDim[0]]), \
        np.array([0,1,0]), np.array([1,0,0])))
        # Step back
        self.rectangles.append(Rectangle(origin+self.origin[1], np.array([stepDim[2],stepDim[0]]), \
        np.array([0,0,-1]), np.array([1,0,0])))
        # Left wall
        self.rectangles.append(Rectangle(origin+self.origin[2], np.array([Dim[2]-stepDim[2],stepDim[1]]), \
        np.array([0,0,1]), np.array([0,1,0])))
        # Right wall
        self.rectangles.append(Rectangle(origin+self.origin[3], np.array([Dim[2]-stepDim[2],stepDim[1]]), \
        np.array([0,0,-1]), np.array([0,1,0])))
        # Top wall
        self.rectangles.append(Rectangle(origin+self.origin[4], np.array([stepDim[0], stepDim[1]]), \
        np.array([1,0,0]), np.array([0,1,0])))

    def write(self,file):
        for rectangle in self.rectangles:
            rectangle.write(file)
        for triangle in self.triangles:
            triangle.write(file)

class WedgeElement():
    """Wedge element"""
    def __init__(self, origin, Dim, stepDim):
        self.rectangles = []
        self.triangles = []
        self.origin = np.array( [[0,0,stepDim[2]],\
                                [0,stepDim[1]/2,stepDim[2]],\
                                [stepDim[0],0,Dim[2]],\
                                [stepDim[0],stepDim[1]/2,Dim[2]],\
                                [0,0,Dim[2]],\
                                [0,stepDim[1]/2,Dim[2]]])

        # The wedge vertices are ordered as follows:
        #       1-------4
        #     /  \     / \
        #   /    2---/---5
        # 0 --------3
        wedge = np.array([[0,0,0],[0,0.5*stepDim[1],stepDim[2]],[0,stepDim[1],0],\
                        [stepDim[0],0,0],[stepDim[0],0.5*stepDim[1],stepDim[2]],[stepDim[0],stepDim[1],0]])
        t1Front = np.array([0,0.5*stepDim[1],stepDim[2]])
        t1Back = np.array([0,0.5*stepDim[1],-stepDim[2]])
        t2 = np.array([1,0,0])
        nFront = np.cross(t1Front,t2)
        nFront_hat = nFront / (nFront**2).sum()**0.5
        nBack = np.cross(t1Back,t2)
        nBack_hat = nBack / (nBack**2).sum()**0.5

        for vertex in wedge:
            vertex += origin

        # Wedge front
        self.triangles.append(Triangle(np.array([wedge[0],wedge[1],wedge[3]]),nFront_hat))
        self.triangles.append(Triangle(np.array([wedge[1],wedge[4],wedge[3]]),nFront_hat))
        # Wedge back
        self.triangles.append(Triangle(np.array([wedge[1],wedge[5],wedge[4]]),nBack_hat))
        self.triangles.append(Triangle(np.array([wedge[1],wedge[2],wedge[5]]),nBack_hat))
        # Left wall
        self.triangles.append(Triangle(np.array([wedge[0],[0,0,stepDim[2]]+origin,wedge[1]]),np.array([-1,0,0])))
        self.triangles.append(Triangle(np.array([wedge[1],[0,stepDim[1],stepDim[2]]+origin,wedge[2]]),np.array([-1,0,0])))
        self.rectangles.append(Rectangle(origin+self.origin[0], np.array([Dim[2]-stepDim[2],stepDim[1]/2]), \
        np.array([0,0,1]), np.array([0,1,0])))
        self.rectangles.append(Rectangle(origin+self.origin[1], np.array([Dim[2]-stepDim[2],stepDim[1]/2]), \
        np.array([0,0,1]), np.array([0,1,0])))
        # Right wall
        self.triangles.append(Triangle(np.array([wedge[3],wedge[4],[stepDim[0],0,stepDim[2]]+origin]),np.array([1,0,0])))
        self.triangles.append(Triangle(np.array([wedge[5],[stepDim[0],stepDim[1],stepDim[2]]+origin,wedge[4]]),np.array([1,0,0])))
        self.rectangles.append(Rectangle(origin+self.origin[2], np.array([Dim[2]-stepDim[2],stepDim[1]/2]), \
        np.array([0,0,-1]), np.array([0,1,0])))
        self.rectangles.append(Rectangle(origin+self.origin[3], np.array([Dim[2]-stepDim[2],stepDim[1]/2]), \
        np.array([0,0,-1]), np.array([0,1,0])))
        # Top wall
        self.rectangles.append(Rectangle(origin+self.origin[4], np.array([stepDim[0], stepDim[1]/2]), \
        np.array([1,0,0]), np.array([0,1,0])))
        self.rectangles.append(Rectangle(origin+self.origin[5], np.array([stepDim[0], stepDim[1]/2]), \
        np.array([1,0,0]), np.array([0,1,0])))

    def write(self,file):
        for rectangle in self.rectangles:
            rectangle.write(file)
        for triangle in self.triangles:
            triangle.write(file)

class RampElement():
    """Ramp element"""
    def __init__(self,origin,Dim,stepDim):
        self.rectangles = []
        self.triangles = []
        self.origin = np.array([[0,stepDim[1],stepDim[2]],\
                                [0,0,stepDim[2]],\
                                [stepDim[0],0,Dim[2]],\
                                [0,0,Dim[2]]])
        # The vertex ordering is the same as for the wedge
        ramp = np.array([[0,0,0],[0,stepDim[1],stepDim[2]],[0,stepDim[1],0],\
                        [stepDim[0],0,0],[stepDim[0],stepDim[1],stepDim[2]],[stepDim[0],stepDim[1],0] ])
        t1 = np.array([0,stepDim[1],stepDim[2]])
        t2 = np.array([1,0,0])
        n = np.cross(t1,t2)
        n_hat = n / (n**2).sum()**0.5

        for vertex in ramp:
            vertex += origin

        # Ramp front (diagonal)
        self.triangles.append(Triangle(np.array([ramp[0],ramp[1],ramp[3]]),n))
        self.triangles.append(Triangle(np.array([ramp[1],ramp[4],ramp[3]]),n))
        # Ramp back (vertical)
        self.rectangles.append(Rectangle(origin+self.origin[0],np.array([stepDim[2],stepDim[0]]), \
        np.array([0,0,-1]),np.array([1,0,0])))
        # Left wall
        self.triangles.append(Triangle(np.array([ramp[0],[0,0,stepDim[2]]+origin,ramp[1]]),np.array([-1,0,0])))
        self.rectangles.append(Rectangle(origin+self.origin[1], np.array([Dim[2]-stepDim[2],stepDim[1]]), \
        np.array([0,0,1]), np.array([0,1,0])))
        # Right wall
        self.triangles.append(Triangle(np.array([ramp[3],ramp[4],[stepDim[0],0,stepDim[2]]+origin]),np.array([1,0,0])))
        self.rectangles.append(Rectangle(origin+self.origin[2], np.array([Dim[2]-stepDim[2],stepDim[1]]), \
        np.array([0,0,-1]), np.array([0,1,0])))
        # Top wall
        self.rectangles.append(Rectangle(origin+self.origin[3], np.array([stepDim[0], stepDim[1]]), \
        np.array([1,0,0]), np.array([0,1,0])))

    def write(self,file):
        for rectangle in self.rectangles:
            rectangle.write(file)
        for triangle in self.triangles:
            triangle.write(file)

class ReverseRampElement():
    """Reverse ramp element"""
    def __init__(self,origin,Dim,stepDim):
        self.rectangles = []
        self.triangles = []
        self.origin = np.array([[0,stepDim[1],stepDim[2]],\
                                [0,0,stepDim[2]],\
                                [stepDim[0],0,Dim[2]],\
                                [0,0,Dim[2]]])
        # The vertex ordering is the same as for the wedge
        ramp = np.array([[0,0,0],[0,0,stepDim[2]],[0,stepDim[1],0],\
                        [stepDim[0],0,0],[stepDim[0],0,stepDim[2]],[stepDim[0],stepDim[1],0] ])
        t1 = np.array([0,stepDim[1],-stepDim[2]])
        t2 = np.array([1,0,0])
        n = np.cross(t1,t2)
        n_hat = n / (n**2).sum()**0.5

        for vertex in ramp:
            vertex += origin

        # Ramp front (vertical)
        self.rectangles.append(Rectangle(origin,np.array([stepDim[2],stepDim[0]]), \
        np.array([0,0,1]),np.array([1,0,0])))
        # Ramp back (diagonal)
        self.triangles.append(Triangle(np.array([ramp[1],ramp[5],ramp[4]]),n))
        self.triangles.append(Triangle(np.array([ramp[1],ramp[2],ramp[5]]),n))
        # Left wall
        self.triangles.append(Triangle(np.array([ramp[2],ramp[1],[0,stepDim[1],stepDim[2]]+origin]),np.array([-1,0,0])))
        self.rectangles.append(Rectangle(origin+self.origin[1], np.array([Dim[2]-stepDim[2],stepDim[1]]), \
        np.array([0,0,1]), np.array([0,1,0])))
        # Right wall
        self.triangles.append(Triangle(np.array([ramp[5],[stepDim[0],stepDim[1],stepDim[2]]+origin,ramp[4]]),np.array([1,0,0])))
        self.rectangles.append(Rectangle(origin+self.origin[2], np.array([Dim[2]-stepDim[2],stepDim[1]]), \
        np.array([0,0,-1]), np.array([0,1,0])))
        # Top wall
        self.rectangles.append(Rectangle(origin+self.origin[3], np.array([stepDim[0], stepDim[1]]), \
        np.array([1,0,0]), np.array([0,1,0])))

    def write(self,file):
        for rectangle in self.rectangles:
            rectangle.write(file)
        for triangle in self.triangles:
            triangle.write(file)

class CubeElement():
    """ A cubic element, used for 3D channels with multiple cubes."""

    # Vertex order, first channel floor, then the cube element
    #          9------10
    #       1-/|     /|----2
    #      / 8------11|   /
    #     /  | 5----|-6  /
    #    /   |/     |/  /
    #   /    4------7  /
    #  0--------------3

    def __init__(self,origin,Dim,cubeOrigin,cubeDim):
        self.vertices = []
        self.triangles = []
        self.rectangles = []

        self.origin = origin
        self.cubeOrigin = cubeOrigin

        self.vertices = np.array([\
                                [0,0,0],\
                                [0,Dim[1],0],\
                                [Dim[0],Dim[1],0],\
                                [Dim[0],0,0],\
                                [cubeOrigin[0],cubeOrigin[1],0],\
                                [cubeOrigin[0],cubeOrigin[1]+cubeDim[1],0],\
                                [cubeOrigin[0]+cubeDim[0],cubeOrigin[1]+cubeDim[1],0],\
                                [cubeOrigin[0]+cubeDim[0],cubeOrigin[1],0],\
                                [cubeOrigin[0],cubeOrigin[1],cubeDim[2]],\
                                [cubeOrigin[0],cubeOrigin[1]+cubeDim[1],cubeDim[2]],\
                                [cubeOrigin[0]+cubeDim[0],cubeOrigin[1]+cubeDim[1],cubeDim[2]],\
                                [cubeOrigin[0]+cubeDim[0],cubeOrigin[1],cubeDim[2]]])

        for vertex in self.vertices:
            vertex += self.origin

        # Cube Front
        self.rectangles.append(Rectangle(self.vertices[4], np.array([cubeDim[2],cubeDim[0]]), np.array([0,0,1]), np.array([1,0,0])))
        # Cube Back
        self.rectangles.append(Rectangle(self.vertices[9], np.array([cubeDim[2],cubeDim[0]]), np.array([0,0,-1]), np.array([1,0,0])))
        # Cube Left
        self.rectangles.append(Rectangle(self.vertices[4], np.array([cubeDim[1],cubeDim[2]]), np.array([0,1,0]), np.array([0,0,1])))
        # Cube Right
        self.rectangles.append(Rectangle(self.vertices[7], np.array([cubeDim[2],cubeDim[1]]), np.array([0,0,1]), np.array([0,1,0])))
        # Cube Top
        self.rectangles.append(Rectangle(self.vertices[8], np.array([cubeDim[1],cubeDim[0]]), np.array([0,1,0]), np.array([1,0,0])))
        # Channel Bottom
        self.triangles.append(Triangle(np.array([self.vertices[0],self.vertices[7],self.vertices[3]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[0],self.vertices[4],self.vertices[7]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[1],self.vertices[4],self.vertices[0]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[1],self.vertices[5],self.vertices[4]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[2],self.vertices[5],self.vertices[1]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[2],self.vertices[6],self.vertices[5]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[3],self.vertices[6],self.vertices[2]]),np.array([0,0,-1])))
        self.triangles.append(Triangle(np.array([self.vertices[3],self.vertices[7],self.vertices[6]]),np.array([0,0,-1])))

    def write(self,file):
        for rectangle in self.rectangles:
            rectangle.write(file)
        for triangle in self.triangles:
            triangle.write(file)

class SubChannel():
    """A rectangular channel object, used to define the regular channels between channel elements."""
    def __init__(self, origin, Dim, stepDim=np.zeros(3)):
        self.rectangles = []
        self.origin = np.array( [[0,0,stepDim[2]],\
                                [0,0,Dim[2]],\
                                [stepDim[0],0,stepDim[2]],\
                                [stepDim[0],0,Dim[2]]])
        # Left wall
        self.rectangles.append(Rectangle(origin, np.array([stepDim[2],Dim[1]]), \
        np.array([0,0,1]), np.array([0,1,0])))
        self.rectangles.append(Rectangle(origin+self.origin[0], np.array([Dim[2]-stepDim[2],Dim[1]]), \
        np.array([0,0,1]), np.array([0,1,0])))
        # Top wall
        self.rectangles.append(Rectangle(origin+self.origin[1], np.array([stepDim[0], Dim[1]]), \
        np.array([1,0,0]), np.array([0,1,0])))
        # Right wall
        self.rectangles.append(Rectangle(origin+self.origin[2], np.array([stepDim[2],Dim[1]]), \
        np.array([0,0,-1]), np.array([0,1,0])))
        self.rectangles.append(Rectangle(origin+self.origin[3], np.array([Dim[2]-stepDim[2],Dim[1]]), \
        np.array([0,0,-1]), np.array([0,1,0])))
        # Bottom wall
        self.rectangles.append(Rectangle(origin, np.array([Dim[1],stepDim[0]]), \
        np.array([0,1,0]), np.array([1,0,0])))

    def write(self,file):
        for rectangle in self.rectangles:
            rectangle.write(file)

class Rectangle():
    """A rectangular object"""

    def __init__(self,origin,dim,t1,t2):
        vertices = np.zeros((4,3))
        vertices[1] = np.array([dim[0]*t1[0], dim[0]*t1[1], dim[0]*t1[2]])
        vertices[3] = np.array([dim[1]*t2[0], dim[1]*t2[1], dim[1]*t2[2]])
        vertices[2] = vertices[1] + vertices[3]
        for vertex in vertices:
            vertex += origin
        n = np.cross(t1,t2)
        # Rectangular vertices
        #   1---2
        #   |   |
        #   0---3
        self.triangle1 = Triangle(np.array([vertices[0],vertices[1],vertices[3]]),n)
        self.triangle2 = Triangle(np.array([vertices[1],vertices[2],vertices[3]]),n)

    def write(self,file):
        self.triangle1.write(file)
        self.triangle2.write(file)

class Triangle():
    """A triangle object"""
    def __init__(self,vertices,n):
        self.vertices = vertices
        self.n = n

    def write(self,file):
        writeString = """facet normal {:.6e} {:.6e} {:.6e}
    outer loop
    vertex {:.6e} {:.6e} {:.6e}
    vertex {:.6e} {:.6e} {:.6e}
    vertex {:.6e} {:.6e} {:.6e}
    endloop
endfacet\n""".format(\
        self.n[0], self.n[1], self.n[2], \
        self.vertices[0][0], self.vertices[0][1], self.vertices[0][2], \
        self.vertices[1][0], self.vertices[1][1], self.vertices[1][2], \
        self.vertices[2][0], self.vertices[2][1], self.vertices[2][2])
        file.write(writeString)
