from Tkinter import *
from math import *
from random import *
from time import *

from scipy.spatial import Delaunay
from numpy import *

class triangulationClass() :
    root = Tk()
    
    trianglationDone = False
    points = []
    faces = []
    tree = []
    
    def __init__(self, _width = 600, _height = 400) :
        self.width = _width
        self.height = _height
        self.root.title("Construction of triangulation")
        self.root.minsize(width = _width, height = _height)
        self.root.maxsize(width = _width, height = _height)

        self.canvas = Canvas(self.root, width = _width, height = _height, bg = "white")       
        self.canvas.pack()

        self.root.bind('<Button-1>', self.clickMouse)

        self.b1 = Button(self.root, bg = "white", fg = "blue", text = "triangulation", command = self.makeTriangle)
        self.b1.place(x = 50, y = _height - 50)        
        
        self.b2 = Button(self.root, bg = "white", fg = "blue", text = "MST", command = self.makeMST)
        self.b2.place(x = 200, y = _height - 50)
        
    def clickMouse(self, event) :
        getX = event.x_root
        getY = event.y_root

        posRootX = self.root.winfo_rootx()
        posRootY = self.root.winfo_rooty()
    
        x = getX - posRootX
        y = getY - posRootY

        self.points.append([x, y])
        self.drawPoint(x, y)
    
    def drawPoint(self, x, y) :
        rad = 1
        if self.trianglationDone == False :
            circle = self.canvas.create_oval(x - rad, y - rad, x + rad, y + rad, fill = "black", outline = "black")
        else :
            self.trianglationDone = False
            self.faces = []
            
            self.canvas.delete("all")
            
            for p in self.points :
                x = p[0]
                y = p[1]
                circle = self.canvas.create_oval(x - rad, y - rad, x + rad, y + rad, fill = "black", outline = "black")
            
    def drawTriangle(self) :
        for i in range(len(self.faces)) :
            temp = [self.points[self.faces[i][0]], self.points[self.faces[i][1]], self.points[self.faces[i][2]], self.points[self.faces[i][0]]]
            for j in range(3) :
                line = self.canvas.create_line(temp[j][0], temp[j][1], temp[j + 1][0], temp[j + 1][1])

    def drawTree(self) :
        for i in range(len(self.tree)) :
            p1 = self.points[self.tree[i][0]]
            p2 = self.points[self.tree[i][1]]
            line = self.canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill = "red")

    def makeTriangle(self) :
        if self.trianglationDone == False :
            self.trianglationDone = True
            
            self.points.pop(len(self.points) - 1)
            tri = Delaunay(self.points)
            
            #print(vars(tri))
            
            for f in tri.vertices :
                self.faces.append([f[0], f[1], f[2]])
        
            self.drawTriangle()

    def getLength(self, p1, p2) :
        return ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2) ** (0.5)

    def makeMST(self) :
        self.makeTriangle()
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

        self.drawTree()