"""These are helper functions used to define channel objects that are written
in STL format"""

import numpy as np
from decimal import Decimal
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import triangle as tr

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

        # Vertices placed in an array
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

            # Vertices placed in an array
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
