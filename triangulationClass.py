from Tkinter import *
from math import *
from random import *
from time import *

from scipy.spatial import Delaunay
import numpy as np

class triangulationClass() :
    root = Tk()
    
    trianglationDone = False
    firstTriangle = True

    points = [[], []]
    faces = [[], []]
    tree = [[], []]

    color = ["#FFCC66", "#0099FF"]
    colorMST = ["#FF6600", "#000099"]
    
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

            self.points[self.firstTriangle == False].append([x, y])
            
            self.drawPoint(x, y)
    
    def secondTriangle(self) :
        self.firstTriangle = False

    def drawPoint(self, x, y) :
        rad = 2
        if self.trianglationDone == False :
            c = self.color[self.firstTriangle == False]
            circle = self.canvas.create_oval(x - rad, y - rad, x + rad, y + rad, fill = c, outline = c)
            
    def drawTriangle(self, num) :
        for i in range(len(self.faces[num])) :
            f = self.faces[num][i]
            temp = [self.points[num][f[0]], self.points[num][f[1]], self.points[num][f[2]], self.points[num][f[0]]]
            
            for j in range(3) :
                line = self.canvas.create_line(temp[j][0], temp[j][1], temp[j + 1][0], temp[j + 1][1], fill = self.color[num])

    def makeTriangle(self) :
        if self.trianglationDone == False :
            self.trianglationDone = True
            
            for i in range(2) :
                self.points[i].pop(len(self.points[i]) - 1)
                
                tri = Delaunay(self.points[i])

                for f in tri.vertices :
                    self.faces[i].append([f[0], f[1], f[2]])

                self.drawTriangle(i)
            
    def drawTree(self) :
        for i in range(2) :
            for j in range(len(self.tree[i])) :
                p1 = self.points[i][self.tree[i][j][0]]
                p2 = self.points[i][self.tree[i][j][1]]
            
                line = self.canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill = self.colorMST[i])

    def getLength(self, p1, p2) :
        return ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2) ** (0.5)

    def makeMST(self) :
        self.makeTriangle()
        edges = [[], []]
        self.tree = [[], []]
        
        for i in range(2) :
            for f in self.faces[i] :
                p = [f[0], f[1], f[2], f[0]]
            
                for j in range(3) :
                    lenEdge = self.getLength(self.points[i][p[j]], self.points[i][p[j + 1]])

                    if edges[i].count((lenEdge, p[j], p[j + 1])) == 0 :
                        edges[i].append((lenEdge, p[j], p[j + 1]))

            edges[i].sort()
            
            treeId = []
            for j in range(len(self.points[i])) :
                treeId.append(j)

            for e in edges[i] :
                l = e[0]
                a = e[1]
                b = e[2]

                if treeId[a] != treeId[b] :
                    self.tree[i].append((a, b))
                    oldId = treeId[b]
                    newId = treeId[a]
                    for j in range(len(treeId)) :
                        if treeId[j] == oldId :
                            treeId[j] = newId

        self.drawTree()

    def experiment(self) :
        self.trianglationDone = False
        self.firstTriangle = True
        self.points = [[], []]
        self.faces = [[], []]
        self.canvas.delete("all")

        nPoints = [10, 10]
        for i in range(2) :
            x = np.random.randint(0, self.width - 20, (nPoints[i], 1)) + 10
            y = np.random.randint(0, self.height - 60, (nPoints[i], 1)) + 5
            self.points[i] = np.concatenate((x, y), axis = 1)

            for j in range(len(self.points[i]) - 1) :
                self.drawPoint(self.points[i][j][0], self.points[i][j][1])

            self.firstTriangle = False

        self.points[0] = self.points[0].tolist()
        self.points[1] = self.points[1].tolist()
        self.makeMST()