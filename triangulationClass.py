from Tkinter import *

from math import *
from random import *
from time import *

from scipy.spatial import Delaunay
import numpy as np
import Queue

EPS = 1e-4

class triangulationClass() :
    root = Tk()
    
    startTriDone = False
    checkTriDone = False
    firstTri = True

    points = [[], [], []]
    nPoints = [[], [], []]
    faces = [[], [], []]
    neighFaces = [[], [], []]
    tree = [[], []]

    color = ["#C0DCC0", "#000080", "#00CC33"]
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
        #print self.points[0]
        #print self.points[1]

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

        return res

    def countSignAngle(self, p1, p2, p3) :
        v1 = [p1[0] - p2[0], p1[1] - p2[1]]
        v2 = [p3[0] - p2[0], p3[1] - p2[1]]

        res = v1[0] * v2[1] - v2[0] * v1[1]
        return res

    def createStruct(self) :
        self.neighbors = [[], [], []]

        for i in range(2) :
            self.nPoints[i] = len(self.points[i])
            temp = [[] for x in range(self.nPoints[i])]
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
                    n = len(temp[tempFace[j]])
                    if n == 0 :
                        angle = self.countSignAngle(self.points[i][tempFace[j - 1]], self.points[i][tempFace[j]], self.points[i][tempFace[j + 1]])
                        if angle < 0 :
                            temp[tempFace[j]] = [tempFace[j - 1], tempFace[j + 1]]
                        else :
                            temp[tempFace[j]] = [tempFace[j + 1], tempFace[j - 1]]
                    else :
                        c1 = temp[tempFace[j]].count(tempFace[j - 1])
                        c2 = temp[tempFace[j]].count(tempFace[j + 1])

                        if temp[tempFace[j]][0] == tempFace[j - 1] and c2 == 0 :
                            temp[tempFace[j]].insert(0, tempFace[j + 1])

                        elif temp[tempFace[j]][0] == tempFace[j + 1] and c1 == 0 :
                            temp[tempFace[j]].insert(0, tempFace[j - 1])

                        elif temp[tempFace[j]][n - 1] == tempFace[j - 1] and c2 == 0 :
                            temp[tempFace[j]].insert(n, tempFace[j + 1])

                        elif temp[tempFace[j]][n - 1] == tempFace[j + 1] and c1 == 0 :
                            temp[tempFace[j]].insert(n, tempFace[j - 1])         

            self.neighbors[i] = temp

        self.nPoints[2] = self.nPoints[0] + self.nPoints[1]
        self.neighbors[2] = [[] for i in range(self.nPoints[2])]

        for i in range(2) :
            for j in range(self.nPoints[i]) :
                for k in range(len(self.neighbors[i][j])) :
                    self.neighbors[2][i * self.nPoints[0] + j].append(self.neighbors[i][j][k] + i * self.nPoints[0])

    def countAngle(self, p1, p2, p3) :
        l1 = self.getLength(p1, p2)
        l2 = self.getLength(p2, p3)
        
        if l1 * l2 == 0 :
            return 0

        v1 = [p2[0] - p1[0], p2[1] - p1[1]]
        v2 = [p2[0] - p3[0], p2[1] - p3[1]]

        aSin = asin((v1[0] * v2[1] - v2[0] * v1[1]) / (l1 * l2))
        aCos = acos((v1[0] * v2[0] + v1[1] * v2[1]) / (l1 * l2))
        if aSin >= 0 :
            return aCos
        else :
            return 2 * pi - aCos

    def addPointToPencil(self, where, p_num) :
        if self.neighbors[2][where].count(p_num) > 0 :
            return

        numPoints = len(self.neighbors[2][where])
        resAngle = []
        for j in range(numPoints) :
            res = self.countAngle(self.points[2][self.neighbors[2][where][j]], self.points[2][where], self.points[2][p_num])
            resAngle.append(res)
        
        aMin = resAngle.index(min(resAngle))
        
        self.neighbors[2][where].insert(aMin, p_num)
        
    def addStarterPointsToPencil(self, starter) :
        for i in range(2) :
            self.addPointToPencil(starter[i], starter[1 - i])

    def deleteWrongEdges(self, starter, reverse = 0) :
        [a_num, b_num] = starter
        a = self.points[2][a_num]
        b = self.points[2][b_num]
        while(1) :
            idx = self.neighbors[2][a_num].index(b_num)
            numPoints = len(self.neighbors[2][a_num])
            c1_num = self.neighbors[2][a_num][(idx - 1) % numPoints]
            c1 = self.points[2][c1_num]
            c2_num = self.neighbors[2][a_num][(idx - 2) % numPoints]
            c2 = self.points[2][c2_num]
            
            """self.drawPoint(a, "purple")
            self.drawPoint(c1, "purple")
            self.drawPoint(c2, "purple")"""

            if c1_num == idx or c2_num == idx:
                break
                
            abc1 = self.countAngle(c1, b, a)
            ac1c2 = self.countAngle(a, c2, c1)
            #print "abc1 ", abc1, " ac1c2 ", ac1c2
            #return 
            
            if abc1 + ac1c2 <= pi or abc1 > pi or ac1c2 > pi:
                break
            else :
                self.neighbors[2][a_num].remove(c1_num)
                self.neighbors[2][c1_num].remove(a_num)
            
        while(1) :
            idx = self.neighbors[2][a_num].index(b_num)
            numPoints = len(self.neighbors[2][a_num])
            c1_num = self.neighbors[2][a_num][(idx + 1) % numPoints]
            c1 = self.points[2][c1_num]
            c2_num = self.neighbors[2][a_num][(idx + 2) % numPoints]
            c2 = self.points[2][c2_num]

            if c1_num == idx or c2_num == idx:
                break
                
            abc1 = self.countAngle(a, b, c1)
            ac1c2 = self.countAngle(c1, c2, a)
            
            if abc1 + ac1c2 <= pi or abc1 > pi or ac1c2 > pi:
                break
            else :
                self.neighbors[2][a_num].remove(c1_num)
                self.neighbors[2][c1_num].remove(a_num)

        if reverse == 0 :
            self.deleteWrongEdges([starter[1], starter[0]], 1)

    def checkNewEdge(self, e1, e2) :
        return (e1[0] == e2[0] and e1[1] == e2[1]) or (e1[0] == e2[1] and e1[1] == e2[0])

    def sewTriangle(self, starter, reverse = 0) :
        [a_num, b_num] = starter
        a = self.points[2][a_num]
        b = self.points[2][b_num]
        
        while(1) :
            idx1 = self.neighbors[2][a_num].index(b_num)
            c_num = self.neighbors[2][a_num][(idx1 - 1) % len(self.neighbors[2][a_num])]
            c = self.points[2][c_num]

            idx2 = self.neighbors[2][b_num].index(a_num)
            d_num = self.neighbors[2][b_num][(idx2 + 1) % len(self.neighbors[2][b_num])]
            d = self.points[2][d_num]
            
            """if a_num == 2 :
                print self.neighbors[2][5]

            print a_num, b_num
            self.drawPoint(a, "purple")
            self.drawPoint(c, "black")
            self.drawPoint(b, "purple")
            self.drawPoint(d, "white")
            
            nb = raw_input()"""

            acb = self.countAngle(a, c, b)
            abd = self.countAngle(a, d, b)
            #print "acb ", acb, " abd ", abd
            #return
            #print "wtf"

            if (acb > abd and  acb < pi) or (abd >= pi and acb < pi) : # CB is edge
                self.drawEdge(c, b, "purple", 2)
                
                self.addPointToPencil(c_num, b_num)
                self.addPointToPencil(b_num, c_num)
                
                aInC = self.neighbors[2][c_num].index(a_num)
                i = (aInC - 1) % len(self.neighbors[2][c_num])
                while(1) :
                    if self.neighbors[2][c_num][i] != b_num :
                        self.neighbors[2][self.neighbors[2][c_num][i]].remove(c_num)
                        self.neighbors[2][c_num].pop(i)
                        i = (i - 1) % len(self.neighbors[2][c_num])
                    else :
                        break

                cInB = self.neighbors[2][b_num].index(c_num)
                i = (cInB - 1) % len(self.neighbors[2][b_num])
                while(1) :
                    if self.neighbors[2][b_num][i] != a_num :
                        self.neighbors[2][self.neighbors[2][b_num][i]].remove(b_num)
                        self.neighbors[2][b_num].pop(i)
                        i = (i - 1) % len(self.neighbors[2][b_num])
                    else :
                        break

                if self.checkNewEdge(starter, [c_num, b_num]) :
                    reverse = 2
                    break
                else :
                    a_num = c_num
                    a = c
            elif (acb < abd and  abd < pi) or (acb >= pi and abd < pi) : # AD is edge
                self.drawEdge(a, d, "yellow", 2)
                
                self.addPointToPencil(d_num, a_num)
                self.addPointToPencil(a_num, d_num)

                aInD = self.neighbors[2][d_num].index(a_num)
                i = (aInD - 1) % len(self.neighbors[2][d_num])
                while(1) :
                    if self.neighbors[2][d_num][i] != b_num :
                        self.neighbors[2][self.neighbors[2][d_num][i]].remove(d_num)
                        self.neighbors[2][d_num].pop(i)
                        i = (i - 1) % len(self.neighbors[2][d_num])
                    else :
                        break

                bInA = self.neighbors[2][a_num].index(b_num)
                i = (bInA - 1) % len(self.neighbors[2][a_num])
                while(1) :
                    if self.neighbors[2][a_num][i] != d_num :
                        self.neighbors[2][self.neighbors[2][a_num][i]].remove(a_num)
                        self.neighbors[2][a_num].pop(i)
                        i = (i - 1) % len(self.neighbors[2][a_num])
                    else :
                        break

                if self.checkNewEdge(starter, [a_num, d_num]) :
                    reverse = 2
                    break
                else :
                    b_num = d_num
                    b = d
            else :
                break
            #nb = raw_input()

        if reverse == 0 :
            self.sewTriangle([starter[1], starter[0]], 1)

    def makeConcatenation(self) :
        self.createStruct()
        
        starter = self.findFirstStarter()
        self.drawEdge(self.points[0][starter[0]], self.points[1][starter[1]], "red")

        starter = [starter[0], starter[1] + self.nPoints[0]]
        self.addStarterPointsToPencil(starter)
        
        self.deleteWrongEdges(starter)
        
        """b_num = starter[0]
        for p in self.neighbors[2][b_num] :
            self.drawPoint(self.points[2][p], "green")"""
        
        """print "light green ", self.neighbors[2][starter[0]]
        print "blue", self.neighbors[2][starter[1]]
        print self.neighbors[2][5]"""
        self.sewTriangle(starter)
        
    def experiment(self) :
        self.startTriDone = False
        self.checkTriDone = False
        self.firstTri = True
        self.points = [[], [], []]
        self.faces = [[], [], []]
        self.canvas.delete("all")

        randomPoints = 0 # 0 - not rand, 1 - rand, 2 - grid rand
        if randomPoints == 0 :
            #self.points[0] = [[147, 93], [226, 45], [253, 130], [149, 197], [78, 54], [-1, -1]]
            #self.points[1] = [[391, 124], [478, 133], [432, 181], [-1, -1]]
            
            self.points[0] = [[80, 209], [155, 141], [223, 183], [217, 264], [138, 297], [92, 365]]
            self.points[1] = [[159, 235], [268, 314], [314, 266], [275, 370]]
            for i in range(2) :
                for j in range(len(self.points[i]) - 1) :
                    self.drawPoint(self.points[i][j], self.color[i])
        elif randomPoints == 1 :
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
        else :
            pass

        self.makeTriangle()