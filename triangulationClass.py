from Tkinter import *
from math import *
from random import *
from time import *

from scipy.spatial import Delaunay
import numpy as np

color1 = "#FFCC66"
color2 = "#0099FF"

class triangulationClass() :
    root = Tk()
    
    trianglationDone = False
    firstTriangle = True

    points1 = []
    points2 = []

    faces1 = []
    faces2 = []

    #tree = []
    
    
    def __init__(self, _width = 600, _height = 400) :
        self.width = _width
        self.height = _height
        self.root.title("Construction of triangulation")
        self.root.minsize(width = _width, height = _height)
        self.root.maxsize(width = _width, height = _height)

        self.canvas = Canvas(self.root, width = _width, height = _height, bg = "white")       
        self.canvas.pack()

        self.root.bind('<Button-1>', self.clickMouse)

        self.b1 = Button(self.root, bg = "white", fg = "blue", text = "Second points", command = self.secondTriangle)
        self.b1.place(x = 50, y = _height - 50)        
        
        self.b2 = Button(self.root, bg = "white", fg = "blue", text = "Triangulation", command = self.makeTriangle)
        self.b2.place(x = 200, y = _height - 50)        
        
        self.b3 = Button(self.root, bg = "white", fg = "blue", text = "MST", command = self.makeMST)
        self.b3.place(x = 350, y = _height - 50)

        self.b4 = Button(self.root, bg = "white", fg = "blue", text = "Experiment", command = self.experiment)
        self.b4.place(x = 450, y = _height - 50)
        
        
    def clickMouse(self, event) :
        if self.trianglationDone == False :
            getX = event.x_root
            getY = event.y_root

            posRootX = self.root.winfo_rootx()
            posRootY = self.root.winfo_rooty()
    
            x = getX - posRootX
            y = getY - posRootY

            if self.firstTriangle == True :
                self.points1.append([x, y])
            else :
                self.points2.append([x, y])
            
            self.drawPoint(x, y)
    
    def secondTriangle(self) :
        self.firstTriangle = False

    def drawPoint(self, x, y) :
        rad = 2
        if self.trianglationDone == False :
            if self.firstTriangle == True :
                circle = self.canvas.create_oval(x - rad, y - rad, x + rad, y + rad, fill = color1, outline = color1)
            else :
                circle = self.canvas.create_oval(x - rad, y - rad, x + rad, y + rad, fill = color2, outline = color2)
        """else :
            self.trianglationDone = False
            self.faces = []
            
            self.canvas.delete("all")
            
            for p in self.points :
                x = p[0]
                y = p[1]
                circle = self.canvas.create_oval(x - rad, y - rad, x + rad, y + rad, fill = "black", outline = "black")"""
            
    def drawTriangle(self, num) :
        if num == 1 :
            for i in range(len(self.faces1)) :
                f = self.faces1[i]
                temp = [self.points1[f[0]], self.points1[f[1]], self.points1[f[2]], self.points1[f[0]]]
                
                for j in range(3) :
                    line = self.canvas.create_line(temp[j][0], temp[j][1], temp[j + 1][0], temp[j + 1][1], fill = color1)
        elif num == 2 :
            for i in range(len(self.faces2)) :
                f = self.faces2[i]
                temp = [self.points2[f[0]], self.points2[f[1]], self.points2[f[2]], self.points2[f[0]]]
                
                for j in range(3) :
                    line = self.canvas.create_line(temp[j][0], temp[j][1], temp[j + 1][0], temp[j + 1][1], fill = color2)

    def makeTriangle(self) :
        self.points1 = np.array(self.points1)
        self.points2 = np.array(self.points2)
        
        if self.trianglationDone == False :
            self.trianglationDone = True
            
            self.points1 = np.delete(self.points1, [len(self.points1) - 1], axis = 0)
            print self.points1
            tri1 = Delaunay(self.points1)

            for f in tri1.vertices :
                self.faces1.append([f[0], f[1], f[2]])

            self.drawTriangle(1)
            

            self.points2 = np.delete(self.points2, [len(self.points2) - 1], axis = 0)
            tri2 = Delaunay(self.points2)

            for f in tri2.vertices :
                self.faces2.append([f[0], f[1], f[2]])            
        
            self.drawTriangle(2)            

    def experiment(self) :
        self.trianglationDone = False
        self.firstTriangle = True
        self.points1 = []
        self.points2 = []
        self.faces1 = []
        self.faces2 = []
        self.canvas.delete("all")

        nPoints = 5
        x = np.random.randint(0, self.width - 20, (nPoints, 1)) + 10
        y = np.random.randint(0, self.height - 60, (nPoints, 1)) + 5
        self.points1 = np.concatenate((x, y), axis = 1)

        for i in range(len(self.points1) - 1) :
            self.drawPoint(self.points1[i][0], self.points1[i][1])

        self.firstTriangle = False

        x = np.random.randint(0, self.width - 10, (nPoints, 1)) + 5
        y = np.random.randint(0, self.height - 60, (nPoints, 1)) + 5
        self.points2 = np.concatenate((x, y), axis = 1)

        for i in range(len(self.points2) - 1) :
            self.drawPoint(self.points2[i][0], self.points2[i][1])

        self.makeTriangle()

    """def drawTree(self) :
        for i in range(len(self.tree)) :
            p1 = self.points[self.tree[i][0]]
            p2 = self.points[self.tree[i][1]]
            line = self.canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill = "red")

    def getLength(self, p1, p2) :
        return ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2) ** (0.5)"""

    def makeMST(self) :
        pass
        """self.makeTriangle()
        edges = []
        for f in self.faces :
            p = [f[0], f[1], f[2], f[0]]
            
            for i in range(3) :
                lenEdge = self.getLength(self.points[p[i]], self.points[p[i + 1]])
                if edges.count((lenEdge, p[i], p[i + 1])) == 0 :
                    edges.append((lenEdge, p[i], p[i + 1]))

        edges.sort()

        treeId = []
        for i in range(len(self.points)) :
            treeId.append(i)

        self.tree = []
        for e in edges :
            l = e[0]
            a = e[1]
            b = e[2]

            if treeId[a] != treeId[b] :
                self.tree.append((a, b))
                oldId = treeId[b]
                newId = treeId[a]
                for i in range(len(treeId)) :
                    if treeId[i] == oldId :
                        treeId[i] = newId

        self.drawTree()"""