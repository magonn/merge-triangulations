from Tkinter import *
from math import *
from random import *
from time import *

from scipy.spatial import Delaunay
import numpy as np

class triangulationClass() :
    root = Tk()
    
    startTriDone = False
    checkTriDone = False
    firstTri = True

    points = [[], [], []]
    faces = [[], [], []]
    tree = [[], []]

    color = ["#FFCC66", "#0099FF", "#00CC33"]
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
        if self.startTriDone == False :
            getX = event.x_root
            getY = event.y_root

            posRootX = self.root.winfo_rootx()
            posRootY = self.root.winfo_rooty()
    
            x = getX - posRootX
            y = getY - posRootY

            p = [x, y]
            self.points[self.firstTri == False].append(p)
            
            self.drawPoint(p, self.color[self.firstTri == False])
    
    def secondTriangle(self) :
        self.firstTri = False

    def drawPoint(self, p, col) :
        rad = 4
        if self.startTriDone == False :
            circle = self.canvas.create_oval(p[0] - rad, p[1] - rad, p[0] + rad, p[1] + rad, fill = col, outline = col)
            
    def drawEdge(self, p1, p2, col) :
        line = self.canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill = col, width = 1)

    def drawTriangle(self, points, faces, col, wid = 1) :
        for f in faces :
            temp = [points[f[0]], points[f[1]], points[f[2]], points[f[0]]]
            
            for j in range(3) :
                line = self.canvas.create_line(temp[j][0], temp[j][1], temp[j + 1][0], temp[j + 1][1], fill = col, width = wid)

    def makeTriangle(self) :
        if self.startTriDone == False :
            self.startTriDone = True
            
            for i in range(2) :
                self.points[i].pop(len(self.points[i]) - 1)
                
                tri = Delaunay(self.points[i])
                
                self.faces[i].extend(tri.vertices)

                self.drawTriangle(self.points[i], self.faces[i], self.color[i])

        if self.checkTriDone == False :
            self.checkTriDone = True

            self.points[2].extend(self.points[0])
            self.points[2].extend(self.points[1])
            
            tri = Delaunay(self.points[2])
            self.faces[2].extend(tri.vertices)

            self.drawTriangle(self.points[2], self.faces[2], self.color[2], 3)

        self.makeConcatenation()

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

    def findFirstStarter(self) :
        starter = [0, 0]
        for i in range(2) :
            for j in range(len(self.points[i])) :
                p = self.points[i][j]
                tempStarter = self.points[i][starter[i]]
                if (p[0] < tempStarter[0]) or (p[0] == tempStarter[0] and p[1] > tempStarter[1]) :
                    starter[i] = j
        
        pStarter = [self.points[0][starter[0]], self.points[1][starter[1]]]
        fixPoint = (pStarter[0][0] < pStarter[1][0]) or (pStarter[0][0] == pStarter[1][0] and (pStarter[0][1] > pStarter[1][1]))

        circle = self.canvas.create_oval(pStarter[fixPoint][0] - 3, pStarter[fixPoint][1] - 3, pStarter[fixPoint][0] + 3, pStarter[fixPoint][1] + 3, fill = "red", outline = "red")
        
        minLength = self.getLength(pStarter[0], pStarter[1])
        for j in range(len(self.points[1 - fixPoint])) :
            length = self.getLength(pStarter[fixPoint], self.points[1 - fixPoint][j])
            if length < minLength :
                minLength = length
                starter[1 - fixPoint] = j
                pStarter[1 - fixPoint] = self.points[1 - fixPoint][j]

        return starter, pStarter

    def makeConcatenation(self) :
        starter, pStarter = self.findFirstStarter()
        self.drawEdge(pStarter[0], pStarter[1], "red")

    def experiment(self) :
        self.startTriDone = False
        self.checkTriDone = False
        self.firstTri = True
        self.points = [[], [], []]
        self.faces = [[], [], []]
        self.canvas.delete("all")

        nPoints = [10, 10]
        for i in range(2) :
            x = np.random.randint(0, self.width - 20, (nPoints[i], 1)) + 10
            y = np.random.randint(0, self.height - 60, (nPoints[i], 1)) + 5
            self.points[i] = np.concatenate((x, y), axis = 1)

            for j in range(len(self.points[i]) - 1) :
                self.drawPoint(self.points[i][j], self.color[i])

            self.secondTriangle()

        self.points[0] = self.points[0].tolist()
        self.points[1] = self.points[1].tolist()
        self.makeTriangle()
        #self.makeMST()
        self.makeConcatenation()