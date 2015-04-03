from Tkinter import *

from math import *
from random import *
from time import *

from scipy.spatial import Delaunay
import numpy as np
import Queue
import copy

class triangulationClass() :
    root = Tk()
    
    startTriDone = False
    checkTriDone = False
    firstTri = True

    points = [[], [], []]
    faces = [[], [], []]
    neighFaces = [[], [], []]
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

    def drawPoint(self, p, col, wid = 4) :
        circle = self.canvas.create_oval(p[0] - wid, p[1] - wid, p[0] + wid, p[1] + wid, fill = col, outline = col)
            
    def drawEdge(self, p1, p2, col, wid = 1) :
        line = self.canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill = col, width = wid)

    def drawTriangle(self, points, faces, col, wid = 1) :
        for f in faces :
            temp = [points[f[0]], points[f[1]], points[f[2]], points[f[0]]]
            
            for j in range(3) :
                line = self.canvas.create_line(temp[j][0], temp[j][1], temp[j + 1][0], temp[j + 1][1], fill = col, width = wid)

    def makeTriangle(self) :
        if self.startTriDone == False :
            self.startTriDone = True
            
            for i in range(2) :
                if len(self.points[i]) < 4 :
                    continue

                self.points[i].pop(len(self.points[i]) - 1)
                
                tri = Delaunay(self.points[i])
                self.faces[i] = tri.vertices
                self.neighFaces[i] = tri.neighbors

                self.drawTriangle(self.points[i], self.faces[i], self.color[i])

        if self.checkTriDone == False :
            self.checkTriDone = True

            self.points[2].extend(self.points[0])
            self.points[2].extend(self.points[1])
            
            tri = Delaunay(self.points[2])
            self.faces[2] = tri.vertices
            self.neighFaces[2] = tri.neighbors

            #self.drawTriangle(self.points[2], self.faces[2], self.color[2], 3)

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
        twoStarter = [0, 0]
        for i in range(2) :
            for j in range(len(self.points[i])) :
                p = self.points[i][j]
                tempStarter = self.points[i][twoStarter[i]]
                if (p[0] < tempStarter[0]) or (p[0] == tempStarter[0] and p[1] > tempStarter[1]) :
                    twoStarter[i] = j
        
        pTwoStarter = [self.points[0][twoStarter[0]], self.points[1][twoStarter[1]]]
        fixPoint = (pTwoStarter[0][0] < pTwoStarter[1][0]) or (pTwoStarter[0][0] == pTwoStarter[1][0] and (pTwoStarter[0][1] > pTwoStarter[1][1]))
        
        endFixPoint = twoStarter[1 - fixPoint]
        pEndFixPoint = pTwoStarter[1 - fixPoint]

        minLength = 1e+10
        for j in range(len(self.points[1 - fixPoint])) :
            p = self.points[1 - fixPoint][j]
            if p[0] == pTwoStarter[fixPoint][0] :
                continue
            
            temp = ((p[0] - pTwoStarter[fixPoint][0]) ** 2 + (p[1] - pTwoStarter[fixPoint][1]) ** 2) / (2 * (pTwoStarter[fixPoint][0] - p[0]))
            if 0 < temp and temp < minLength :
                minLength = temp
                endFixPoint = j
                pEndFixPoint = p
                
        pRes = [pTwoStarter[fixPoint], pEndFixPoint]
        if fixPoint == 1 :
            res = [endFixPoint, twoStarter[1]]
        else :
            res = [twoStarter[0], endFixPoint]

        return res, pRes

    def countSignAngle(self, p1, p2, p3) :
        v1 = [p2[0] - p1[0], p2[1] - p1[1]]
        v2 = [p2[0] - p3[0], p2[1] - p3[1]]

        return v1[0] * v2[1] - v2[0] * v1[1]


    def createStruct(self) :
        self.neighbors = [[], [], []]

        for i in range(2) :
            temp = [[] for x in range(len(self.points[i]))]
            flag = [False for x in range(len(self.faces[i]))]
            
            queue = Queue.Queue()
            queue.put(0)

            while not queue.empty() :
                cur = queue.get()
                
                if flag[cur] == True :
                    continue

                flag[cur] = True
                
                for j in range(3) :
                    if self.neighFaces[i][cur][j] != -1 :
                        queue.put(self.neighFaces[i][cur][j])

                x = self.faces[i][cur]
                tempFace = [x[0], x[1], x[2], x[0], x[1]]
                
                for j in range(1, 4) :
                    if len(temp[tempFace[j]]) == 0 :
                        if self.countSignAngle(self.points[i][tempFace[j - 1]], self.points[i][tempFace[j]], self.points[i][tempFace[j + 1]]) < 0 :
                            temp[tempFace[j]] = [tempFace[j - 1], tempFace[j + 1]]
                        else :
                            temp[tempFace[j]] = [tempFace[j + 1], tempFace[j - 1]]
                    elif temp[tempFace[j]][0] == tempFace[j - 1] and temp[tempFace[j]].count(tempFace[j + 1]) == 0:
                        temp[tempFace[j]].insert(0, tempFace[j + 1])
                    elif temp[tempFace[j]][0] == tempFace[j + 1] and temp[tempFace[j]].count(tempFace[j - 1]) == 0:
                        temp[tempFace[j]].insert(0, tempFace[j - 1])
                    elif temp[tempFace[j]][len(temp[tempFace[j]]) - 1] == tempFace[j - 1] and temp[tempFace[j]].count(tempFace[j + 1]) == 0:
                        temp[tempFace[j]].insert(len(temp[tempFace[j]]), tempFace[j + 1])
                    elif temp[tempFace[j]][len(temp[tempFace[j]]) - 1] == tempFace[j + 1] and temp[tempFace[j]].count(tempFace[j - 1]) == 0 :
                        temp[tempFace[j]].insert(len(temp[tempFace[j]]), tempFace[j - 1])

            self.neighbors[i] = temp
            
        self.neighbors[2] = copy.deepcopy(self.neighbors[0])
        n1 = len(self.neighbors[0])
        for i in range(len(self.points[1])) :
            temp = []
            for j in range(len(self.neighbors[1][i])) :
                temp.append(self.neighbors[1][i][j] + n1)
            self.neighbors[2].append(temp)

    def makeConcatenation(self) :
        starter, pStarter = self.findFirstStarter()
        self.drawEdge(self.points[0][starter[0]], self.points[1][starter[1]], "red")

        self.createStruct()

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