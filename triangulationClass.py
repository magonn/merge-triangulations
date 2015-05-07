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
    experimentMode = False
    MSTmode = False

    points = [[], [], []]
    nPoints = [[], [], []]
    faces = [[], [], []]
    neighFaces = [[], [], []]

    tree = [[], []]
    mst = []
    bridges = []

    #color = ["#FA8072", "#87CEEB", "#00CC33"]
    color = ["black", "black", "black"]
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
        
        self.b3 = Button(self.root, bg = "white", fg = "blue", text = "MST", command = self.MST)
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
            
            self.drawPoint(p, self.color[self.firstTri == False], 3, self.firstTri == False)
    
    def secondTriangle(self) :
        self.firstTri = False

    def drawPoint(self, p, col, wid = 2, type = 0) :
        if type == 0 :
            circle = self.canvas.create_oval(p[0] - wid, p[1] - wid, p[0] + wid, p[1] + wid, fill = col, outline = col)
        else :
            circle = self.canvas.create_oval(p[0] - wid, p[1] - wid, p[0] + wid, p[1] + wid, fill = "white", outline = col)

    def drawCircle(self, p, col, wid = 0, type = 0) :
        if type == 0 :
            circle = self.canvas.create_oval(p[0] - wid, p[1] - wid, p[0] + wid, p[1] + wid, outline = col)
        else :
            circle = self.canvas.create_oval(p[0] - wid, p[1] - wid, p[0] + wid, p[1] + wid, outline = col, dash = '- -')

    def drawEdge(self, p1, p2, col, wid = 1, type = 0) :
        if type == 0 or type == 2 :
            line = self.canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill = col, width = wid)
        else :
            line = self.canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill = col, width = wid, dash = '- -')
            
    def drawTriangle(self, points, faces, col, wid = 1, type = 0) :
        for f in faces :
            temp = [points[f[0]], points[f[1]], points[f[2]], points[f[0]]]
            
            for j in range(3) :
                if type == 0 or type == 2:
                    line = self.canvas.create_line(temp[j][0], temp[j][1], temp[j + 1][0], temp[j + 1][1], fill = col, width = wid)
                else :
                    line = self.canvas.create_line(temp[j][0], temp[j][1], temp[j + 1][0], temp[j + 1][1], fill = col, width = wid, dash = '- -')
     
    def makeTriangle(self) :
        sciPyTime = 0
        if self.startTriDone == False :
            self.startTriDone = True
            
            for i in range(2) :
                if self.experimentMode == False :
                    self.points[i].pop(len(self.points[i]) - 1)
                
                tri = Delaunay(self.points[i])
                self.faces[i] = tri.vertices
                self.neighFaces[i] = tri.neighbors

                if self.MSTmode == True :
                    self.drawTriangle(self.points[i], self.faces[i], self.color[i], 1, i)

        if self.checkTriDone == False :
            self.checkTriDone = True
            
            start1 = time()

            self.points[2] = []
            self.points[2].extend(self.points[0])
            self.points[2].extend(self.points[1])
            
            tri = Delaunay(self.points[2])
            self.faces[2] = tri.vertices
            self.neighFaces[2] = tri.neighbors

            finish1 = time()
            sciPyTime = finish1 - start1

            if self.MSTmode == False :
                self.drawTriangle(self.points[2], self.faces[2], "green", 4, 2)
                #self.drawTriangle(self.points[2], self.faces[2], "black", 1, 0)

        if self.MSTmode == False :
            myTime = self.makeConcatenation()
            return [sciPyTime, myTime]
        
    def drawMST(self, num, col = "black", wid = 1) :
        for j in range(len(self.tree[num])) :
            p1 = self.points[num][self.tree[num][j][0]]
            p2 = self.points[num][self.tree[num][j][1]]
            
            line = self.canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill = col, width = wid)

    def getLength(self, p1, p2) :
        return ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2) ** (0.5)

    def makeMST(self) :
        edgesMST = [[], []]
        self.tree = [[], []]
        
        for i in range(2) :
            for f in self.faces[i] :
                p = [f[0], f[1], f[2], f[0]]
            
                for j in range(3) :
                    lenEdge = self.getLength(self.points[i][p[j]], self.points[i][p[j + 1]])

                    if edgesMST[i].count((lenEdge, p[j], p[j + 1])) == 0 :
                        edgesMST[i].append((lenEdge, p[j], p[j + 1]))

            edgesMST[i].sort()
            
            treeId = []
            for j in range(len(self.points[i])) :
                treeId.append(j)

            for e in edgesMST[i] :
                a = e[1]
                b = e[2]

                if treeId[a] != treeId[b] :
                    self.tree[i].append((a, b))
                    oldId = treeId[b]
                    newId = treeId[a]
                    for j in range(len(treeId)) :
                        if treeId[j] == oldId :
                            treeId[j] = newId

            #self.drawMST(i, "black", 2)

        self.mst = []
        n1 = len(self.points[0])
        
        for i in range(2) :
            for j in range(len(self.tree[i])) :
                self.mst.append([self.tree[i][j][0] + i * n1, self.tree[i][j][1] + i * n1])

    def MST(self) :
        self.MSTmode = True
        self.makeTriangle()
        
        self.makeMST()
        self.drawMST(0, "black", 2)
        self.drawMST(1, "black", 2)
        self.drawAllPoints()
        self.MSTmode = False      

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
        
        #draw Empty W circle
        #p0 = pTwoStarter[fixPoint]
        #p1 = pTwoStarter[1 - fixPoint]
        #rad = ((p0[0] - p1[0]) ** 2 + (p0[1] - p1[1]) ** 2) / (2 * (p0[0] - p1[0]))
        #self.drawCircle([p0[0] - rad, p0[1]], "black", rad, 0)

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
        res = [res[0], res[1] + len(self.points[0])]
        return res

    def countSignAngle(self, p1, p2, p3) :
        v1 = [p1[0] - p2[0], p1[1] - p2[1]]
        v2 = [p3[0] - p2[0], p3[1] - p2[1]]

        res = v1[0] * v2[1] - v2[0] * v1[1]
        return res

    def drawStruct(self) :
        for i in range(len(self.points[2])) :
            p1 = self.points[2][i]
            for j in range(len(self.neighbors[2][i])) :
                self.drawEdge(p1, self.points[2][self.neighbors[2][i][j]], "purple", 1, 0)

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
                    if self.neighFaces[i][cur][j] != -1 and flag[self.neighFaces[i][cur][j]] == False:
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
                            c2 = 1

                        elif temp[tempFace[j]][0] == tempFace[j + 1] and c1 == 0 :
                            temp[tempFace[j]].insert(0, tempFace[j - 1])
                            c1 = 1

                        elif temp[tempFace[j]][n - 1] == tempFace[j - 1] and c2 == 0 :
                            temp[tempFace[j]].insert(n, tempFace[j + 1])
                            c2 = 1

                        elif temp[tempFace[j]][n - 1] == tempFace[j + 1] and c1 == 0 :
                            temp[tempFace[j]].insert(n, tempFace[j - 1])
                            c1 = 1

                        if c1 == 0 or c2 == 0 :
                            flag[cur] = False
                            queue.put(cur)

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

        sin_arg = (v1[0] * v2[1] - v2[0] * v1[1]) / (l1 * l2)
        if abs(abs(sin_arg) - 1) < EPS :
            x = copysign(1, sin_arg)
            aSin = asin(x)
        else :
            aSin = asin(sin_arg)

        cos_arg = (v1[0] * v2[0] + v1[1] * v2[1]) / (l1 * l2)
        if abs(abs(cos_arg) - 1) < EPS :
            x = copysign(1, cos_arg)
            aCos = acos(x)
        else :
            aCos = acos(cos_arg)
        
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
        
        if len(resAngle) == 0 :
            aMin = 0
        else :
            aMin = resAngle.index(min(resAngle))
        
        self.neighbors[2][where].insert(aMin, p_num)
        
    def addEdgeToPencil(self, edge) :
        for i in range(2) :
            self.addPointToPencil(edge[i], edge[1 - i])

    def checkBridge(self, edge) :
        if self.mst.count(edge) > 0 or self.mst.count([edge[1], edge[0]]) > 0 :
            self.bridges.append(edge)

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
            
            if c1_num == b_num or c2_num == b_num :
                break
                
            """self.drawPoint(a, "#FAEBD7", 3, 0)
            self.drawPoint(c1, "purple", 3, 0)
            self.drawPoint(c2, "purple", 3, 0)"""
            
            abc1 = self.countAngle(c1, b, a)
            ac2c1 = self.countAngle(a, c2, c1)
            
            """print "a, b, c1, c2", a_num, b_num, c1_num, c2_num
            print "abc1 ", abc1, " ac2c1 ", ac2c1
            nb = raw_input()
            self.drawAllPoints()"""

            if abc1 + ac2c1 <= pi or abc1 > pi or ac2c1 > pi:
                break
            else :
                self.checkBridge([a_num, c1_num])
                self.neighbors[2][a_num].remove(c1_num)
                self.neighbors[2][c1_num].remove(a_num)
                
        while(1) :
            idx = self.neighbors[2][a_num].index(b_num)
            numPoints = len(self.neighbors[2][a_num])
            c1_num = self.neighbors[2][a_num][(idx + 1) % numPoints]
            c1 = self.points[2][c1_num]
            c2_num = self.neighbors[2][a_num][(idx + 2) % numPoints]
            c2 = self.points[2][c2_num]

            if c1_num == b_num or c2_num == b_num :
                break
                
            """self.drawPoint(a, "#FF4500", 3, 0)
            self.drawPoint(c1, "purple", 3, 0)
            self.drawPoint(c2, "purple", 3, 0)"""
            
            abc1 = self.countAngle(a, b, c1)
            ac2c1 = self.countAngle(c1, c2, a)
            
            """print "one a, b, c1, c2", a_num, b_num, c1_num, c2_num
            print "abc1 ", abc1, " ac2c1 ", ac2c1
            nb = raw_input()
            self.drawAllPoints() """

            if abc1 + ac2c1 <= pi or abc1 > pi or ac2c1 > pi:
                break
            else :
                self.checkBridge([a_num, c1_num])
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
        self.drawEdge(a, b, "red", 3)
        
        while(1) :
            self.deleteWrongEdges([a_num, b_num])
            
            idx1 = self.neighbors[2][a_num].index(b_num)
            c_num = self.neighbors[2][a_num][(idx1 - 1) % len(self.neighbors[2][a_num])]
            c = self.points[2][c_num]

            idx2 = self.neighbors[2][b_num].index(a_num)
            d_num = self.neighbors[2][b_num][(idx2 + 1) % len(self.neighbors[2][b_num])]
            d = self.points[2][d_num]
            
            #print "EDGE", a_num, b_num, c_num, d_num
            """self.drawPoint(c, "yellow", 3, 0)
            self.drawPoint(d, "white", 3, 0)
            nb = raw_input()"""

            acb = self.countAngle(a, c, b)
            adb = self.countAngle(a, d, b)

            #print "acb ", acb, "adb ", adb
            
            if (c_num != b_num) and ((acb > adb and acb < pi) or (adb >= pi and acb < pi)) : # CB is edge
                self.addEdgeToPencil([c_num, b_num])
                
                aInC = self.neighbors[2][c_num].index(a_num)
                i = (aInC - 1) % len(self.neighbors[2][c_num])
                while(1) :
                    check = self.neighbors[2][c_num][i]
                    if check != b_num :
                        self.checkBridge([c_num, check])
                        
                        self.neighbors[2][check].remove(c_num)
                        self.neighbors[2][c_num].pop(i)
                        i = (i - 1) % len(self.neighbors[2][c_num])
                    else :
                        break

                cInB = self.neighbors[2][b_num].index(c_num)
                i = (cInB - 1) % len(self.neighbors[2][b_num])
                while(1) :
                    check = self.neighbors[2][b_num][i]
                    if check != a_num :
                        self.checkBridge([b_num, check])

                        self.neighbors[2][check].remove(b_num)
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
            elif (a_num != d_num) and (acb < adb and  adb < pi) or (acb >= pi and adb < pi) : # AD is edge
                self.addEdgeToPencil([a_num, d_num])

                aInD = self.neighbors[2][d_num].index(a_num)
                i = (aInD - 1) % len(self.neighbors[2][d_num])
                while(1) :
                    check = self.neighbors[2][d_num][i]
                    if check != b_num :
                        self.checkBridge([d_num, check])

                        self.neighbors[2][check].remove(d_num)
                        self.neighbors[2][d_num].pop(i)
                        i = (i - 1) % len(self.neighbors[2][d_num])
                    else :
                        break

                bInA = self.neighbors[2][a_num].index(b_num)
                i = (bInA - 1) % len(self.neighbors[2][a_num])
                while(1) :
                    check = self.neighbors[2][a_num][i]
                    if check != d_num :
                        self.checkBridge([a_num, check])

                        self.neighbors[2][check].remove(a_num)
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
            #self.drawEdge(a, b, "purple", 1, 0)
            #nb = raw_input()

        if reverse == 0 :
            self.sewTriangle([starter[1], starter[0]], 1)

    def lineParameters(self, p1, p2) :
        a = p1[1] - p2[1]
        b = p2[0] - p1[0]
        c = p1[0] * p2[1] - p2[0] * p1[1]
        return [a, b, c]

    def perpendicular(self, p, line) :
        [a, b, c] = line
        pa = b
        pb = -a
        pc = -b * p[0] + a * p[1]
        return [pa, pb, pc]

    def intersect(self, line1, line2) :
        [a1, b1, c1] = line1
        [a2, b2, c2] = line2
        det = a1 * b2 - a2 * b1
        if abs(det) < EPS :
            return [1000, 1000]
        x0 = (c2 * b1 - c1 * b2) / det
        y0 = (c2 * a1 - c1 * a2) / (-det)
        return [x0, y0]

    def findNextStarter(self) :
        open = -1
        curBridge = []
        while(open == -1) :
            if len(self.bridges) == 0 :
                return -1

            bridge = self.bridges.pop(0)
            #self.drawEdge(self.points[2][bridge[0]], self.points[2][bridge[1]], "purple", 3, 0)

            for i in range(2) :
                num_p = bridge[i]
                pencil_p = (num_p >= len(self.points[0]))
                flag = True
                for j in range(len(self.neighbors[2][num_p])) :
                    pencil_temp = (self.neighbors[2][num_p][j] >= len(self.points[0]))
                    if pencil_p != pencil_temp :
                        flag = False
                        break
                if flag :
                    curBridge = bridge
                    open = num_p
                    break

        #self.drawPoint(self.points[2][open], "black", 5, 1)

        openPoint = self.points[2][open]
        pencil_p = (open >= len(self.points[0]))
        if pencil_p == 0 :
            start = len(self.points[0])
            end = len(self.points[2])
        else :
            start = 0
            end = len(self.points[0])

        bridgeParameters = self.lineParameters(self.points[2][curBridge[0]], self.points[2][curBridge[1]])

        minLength = 1e+10
        num_endStarter = -1;
        for i in range(start, end) :
            tempPoint = self.points[2][i]
            
            [a, b, c] = self.lineParameters(tempPoint, openPoint)

            midlePoint = [round((tempPoint[0] + openPoint[0]) / 2), round((tempPoint[1] + openPoint[1]) / 2)]
            [pa, pb, pc] = self.perpendicular(midlePoint, [a, b, c])

            p0 = self.intersect(bridgeParameters, [pa, pb, pc])
            """print p0
            self.drawPoint(p0, "black", 5, 0)
            nb = raw_input()"""

            temp = self.getLength(p0, openPoint)
            if min(self.points[2][curBridge[0]][0], self.points[2][curBridge[1]][0]) <= p0[0] and p0[0] <= max(self.points[2][curBridge[0]][0], self.points[2][curBridge[1]][0]) and temp < minLength :
                minLength = temp
                num_endStarter = i

        return [open, num_endStarter]

    def checkStruct(self) :
        flag = False
        for i in range(len(self.points[2])) :
            for j in range(len(self.neighbors[2][i])) :
                num_p = self.neighbors[2][i][j]
                if self.neighbors[2][num_p].count(i) != 1 :
                    flag = True
                    print "bad struct"
                    print i, self.neighbors[2][i]
                    print num_p, self.neighbors[2][num_p]
                    break

        if flag :
            print "There are errors"
            print self.points[0]
            print self.points[1]

    def drawAllPoints(self) :
        for i in range(2) :
            for j in range(len(self.points[i])) :
                self.drawPoint(self.points[i][j], self.color[i], 3, i)

    def makeConcatenation(self) :
        print "self.points[0] =", self.points[0]
        print "self.points[1] =", self.points[1]
        
        self.createStruct()
        self.makeMST()
        
        start = time()        

        starter = self.findFirstStarter()
        self.addEdgeToPencil(starter)
        
        self.sewTriangle(starter)
        #print "bridges", self.bridges
        
        while(len(self.bridges) != 0) :
            starter = self.findNextStarter()
            if starter == -1 :
                break
            self.addEdgeToPencil(starter)
            self.sewTriangle(starter)
        
        finish = time()

        self.drawStruct()
        self.drawAllPoints()
        print "ok"
        return finish - start
        
    def experiment(self) :
        self.startTriDone = False
        self.checkTriDone = False
        self.firstTri = True
        self.points = [[], [], []]
        self.faces = [[], [], []]
        self.bridges = []

        self.canvas.delete("all")
        self.experimentMode = True

        randomPoints = 1 # 0 - example, 1 - rand

        if randomPoints == 0 :
            # separeted triangle seam
            #self.points[0] = [[147, 93], [226, 45], [253, 130], [149, 197], [78, 54]]
            #self.points[1] = [[391, 124], [478, 133], [432, 181]]
            
            #example starter in report
            #self.points[0] = [[80, 209], [155, 141], [223, 183], [217, 264], [138, 297]]
            #self.points[1] = [[159, 235], [268, 314], [314, 266]]
            
            #model example in report
            self.points[0] = [[139, 243], [173, 169], [237, 195], [211, 230], [275, 231], [224, 273], [178, 279], [228, 149], [286, 195], [291, 143], [333, 169], [312, 187], [355, 218], [327, 254], [293, 279], [369, 261], [394, 191], [356, 139]]
            self.points[1] = [[187, 255], [198, 193], [248, 251], [264, 173], [315, 226], [329, 294], [262, 302], [404, 236], [358, 186], [323, 119], [411, 151], [398, 294]]

            # big example
            #self.points[0] = [[144, 161], [272, 247], [70, 341], [360, 205], [474, 294], [87, 153], [351, 127], [543, 136], [468, 324], [447, 38], [571, 129], [176, 182], [179, 341], [452, 212], [95, 320], [31, 236], [265, 75], [548, 21], [364, 191], [111, 78], [369, 88], [454, 105], [42, 233], [49, 125], [503, 175], [56, 93], [92, 248], [474, 92], [19, 269], [467, 283], [427, 239], [463, 15], [232, 20], [405, 141], [135, 177], [483, 37], [440, 102], [184, 160], [553, 147], [338, 330], [121, 26], [306, 152], [137, 93], [107, 124], [448, 339], [481, 77], [397, 113], [570, 285], [11, 38], [343, 113], [379, 39], [16, 146], [238, 196], [481, 125], [24, 207], [97, 81], [535, 194], [114, 46], [80, 213], [512, 268], [72, 276], [185, 21], [429, 278], [90, 212], [75, 65], [512, 7], [287, 5], [532, 284], [371, 234], [455, 311], [367, 320], [31, 151], [406, 188], [359, 252], [457, 334], [479, 185], [39, 213], [236, 298], [217, 29], [123, 275], [508, 285], [294, 67], [375, 219], [462, 118], [196, 117], [15, 179], [318, 9], [525, 290], [121, 86], [47, 207], [321, 27], [539, 206], [74, 82], [33, 43], [360, 144], [50, 101], [392, 160], [550, 126], [38, 33], [159, 110]]
            #elf.points[1] = [[444, 249], [337, 338], [574, 206], [212, 173], [437, 275], [477, 118], [172, 163], [129, 232], [481, 143], [551, 65], [543, 12], [93, 88], [556, 85], [468, 236], [477, 86], [438, 159], [123, 288], [562, 14], [19, 57], [176, 256], [98, 302], [228, 343], [451, 52], [571, 211], [108, 62], [484, 243], [496, 250], [315, 342], [358, 264], [499, 39], [321, 283], [269, 228], [81, 81], [113, 245], [201, 296], [371, 250], [579, 263], [315, 105], [344, 48], [290, 57], [499, 251], [81, 238], [509, 109], [510, 246], [325, 128], [118, 302], [357, 16], [78, 217], [69, 133], [195, 343], [545, 232], [236, 231], [467, 341], [204, 329], [545, 301], [89, 166], [477, 159], [483, 49], [35, 314], [343, 327], [312, 188], [107, 244], [24, 120], [126, 165], [475, 120], [15, 72], [344, 97], [221, 44], [111, 246], [16, 182], [181, 71], [305, 267], [111, 10], [94, 187], [144, 308], [190, 105], [366, 277], [474, 231], [308, 111], [65, 17], [477, 343], [475, 81], [574, 39], [581, 36], [538, 330], [506, 167], [317, 123], [243, 213], [54, 181], [363, 113], [28, 214], [211, 337], [375, 146], [481, 279], [552, 193], [456, 270], [519, 183], [133, 20], [250, 15], [523, 95]]
        elif randomPoints == 1 :
            nPoints = [500, 500]
            for i in range(2) :
                x = np.random.randint(0, self.width - 20, (nPoints[i], 1)) + 10
                y = np.random.randint(0, self.height - 60, (nPoints[i], 1)) + 5
                self.points[i] = np.concatenate((x, y), axis = 1)

                self.secondTriangle()

            self.points[0] = self.points[0].tolist()
            self.points[1] = self.points[1].tolist()
        else :
            p1 = [400, 200]
            p2 = [300, 250]
            self.drawEdge([200, 200], [400, 200], "black", 1, 0)
            self.drawEdge([300, 250], [400, 200], "black", 1, 1)
            self.drawPoint(p1, "black", 3, 0)
            self.drawPoint(p2, "black", 3, 1)
            x = 320
            y = 2 * x - 475
            #line = self.lineParameters(p1, p2)
            #lineP = self.perpendicular([(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2], line)
            self.drawEdge([x, y], [(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2], "black", 1, 1)
            
            self.drawPoint([(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2], "black", 3, 1)
            self.drawPoint([337, 200], "black", 3, 1)
            self.drawPoint([250, 150], "black", 3, 1)
            self.drawPoint([230, 220], "black", 3, 1)
            self.drawPoint([380, 100], "black", 3, 1)
            self.drawCircle([337, 200], "black", 400 - 337, 0)
            
            return

        self.drawAllPoints()
        [sciPyTime, myTime] = self.makeTriangle()
        print "sciPy", sciPyTime, "my", myTime