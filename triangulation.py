from Tkinter import *

from random import *
from time import *

from scipy.spatial import Delaunay
import numpy as np
import Queue

import geo
import draw

class triangulation() :
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

    testMode = False

    #color = ["#FA8072", "#87CEEB", "#00CC33"]
    color = ["black", "black", "black"]
    colorMST = ["#FF6600", "#000099"]
    
    def __init__(self, _width = 600, _height = 400) :
        self.root = Tk()
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
            draw.drawPoint(self.canvas, p, self.color[self.firstTri == False], 3, self.firstTri == False)
    
    def secondTriangle(self) :
        self.firstTri = False

    def makeTriangle(self) :
        #print self.points[0]
        #print self.points[1]

        sciPyTime = 0
        if self.startTriDone == False :
            self.startTriDone = True
            
            for i in range(2) :
                if self.experimentMode == False :
                    self.points[i].pop(len(self.points[i]) - 1)
                
                tri = Delaunay(self.points[i])
                self.faces[i] = tri.vertices
                self.neighFaces[i] = tri.neighbors

                #self.drawTriangle(self.points[i], self.faces[i], self.color[i], 1, i)
                #self.drawAllPoints()

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
                draw.drawTriangle(self.canvas, self.points[2], self.faces[2], "green", 5, 2)
                self.drawAllPoints()
                #self.drawTriangle(self.points[2], self.faces[2], "purple", 2, 0)

        if self.MSTmode == False :
            #myTime = self.makeConcatenation("yellow", 3)
            myTime = self.makeConcatenation()
            #self.testMode = True
            #myTimeTest = self.makeConcatenation()
            return [sciPyTime, myTime]#, myTimeTest]
            
    def drawMST(self, num, col = "black", wid = 1) :
        for j in range(len(self.tree[num])) :
            p1 = self.points[num][self.tree[num][j][0]]
            p2 = self.points[num][self.tree[num][j][1]]
            
            draw.drawEdge(self.canvas, p1, p2, col, wid)

    def makeMST(self) :
        edgesMST = [[], []]
        self.tree = [[], []]
        
        for i in range(2) :
            for f in self.faces[i] :
                p = [f[0], f[1], f[2], f[0]]
            
                for j in range(3) :
                    lenEdge = geo.getLength(self.points[i][p[j]], self.points[i][p[j + 1]])

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

    def drawStruct(self, col = "black", wid = 1) :
        for i in range(len(self.points[2])) :
            p1 = self.points[2][i]
            for j in self.neighbors[2][i] :
                draw.drawEdge(self.canvas, p1, self.points[2][j], col, wid, 0)

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
                        angle = geo.exteriorProd(self.points[i][tempFace[j - 1]], self.points[i][tempFace[j]], self.points[i][tempFace[j + 1]])
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

    def addPointToPencil(self, where, p_num) :
        if self.neighbors[2][where].count(p_num) > 0 :
            return

        numPoints = len(self.neighbors[2][where])
        resAngle = []
        for j in range(numPoints) :
            res = geo.countAngle(self.points[2][self.neighbors[2][where][j]], self.points[2][where], self.points[2][p_num])
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
                
            #self.drawPoint(a, "#FAEBD7", 3, 0)
            #self.drawPoint(c1, "purple", 3, 0)
            #self.drawPoint(c2, "purple", 3, 0)
            
            abc1 = geo.countAngle(c1, b, a)
            ac2c1 = geo.countAngle(a, c2, c1)
            
            #print "a, b, c1, c2", a_num, b_num, c1_num, c2_num
            #print "abc1 ", abc1, " ac2c1 ", ac2c1
            #nb = raw_input()
            #self.drawAllPoints()

            if abc1 + ac2c1 <= geo.PI or abc1 > geo.PI or ac2c1 > geo.PI:
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
                
            #self.drawPoint(a, "#FF4500", 3, 0)
            #self.drawPoint(c1, "purple", 3, 0)
            #self.drawPoint(c2, "purple", 3, 0)
            
            abc1 = geo.countAngle(a, b, c1)
            ac2c1 = geo.countAngle(c1, c2, a)
            
            #print "one a, b, c1, c2", a_num, b_num, c1_num, c2_num
            #print "abc1 ", abc1, " ac2c1 ", ac2c1
            #nb = raw_input()
            #self.drawAllPoints()

            if abc1 + ac2c1 <= geo.PI or abc1 > geo.PI or ac2c1 > geo.PI :
                break
            else :
                self.checkBridge([a_num, c1_num])
                self.neighbors[2][a_num].remove(c1_num)
                self.neighbors[2][c1_num].remove(a_num)
                
        if reverse == 0 :
            self.deleteWrongEdges([starter[1], starter[0]], 1)

    def checkNewEdge(self, e1, e2) :
        return (e1[0] == e2[0] and e1[1] == e2[1]) or (e1[0] == e2[1] and e1[1] == e2[0])

    def sewTriangle(self, starter, reverse = 0, colorEdge = "purple") :
        [b_num, a_num] = starter
        a = self.points[2][a_num]
        b = self.points[2][b_num]
        #self.drawEdge(a, b, colorEdge, 2)
        #self.drawEdge(a, b, "black", 3)
        #self.drawAllPoints()
        #nb = raw_input()
        
        while(1) :
            self.deleteWrongEdges([a_num, b_num])
            
            idx1 = self.neighbors[2][a_num].index(b_num)
            c_num = self.neighbors[2][a_num][(idx1 - 1) % len(self.neighbors[2][a_num])]
            c = self.points[2][c_num]

            idx2 = self.neighbors[2][b_num].index(a_num)
            d_num = self.neighbors[2][b_num][(idx2 + 1) % len(self.neighbors[2][b_num])]
            d = self.points[2][d_num]
            
            #print "EDGE", a_num, b_num, c_num, d_num
            #self.drawPoint(c, "yellow", 3, 0)
            #self.drawPoint(d, "white", 3, 0)
            #nb = raw_input()

            acb = geo.countAngle(a, c, b)
            adb = geo.countAngle(a, d, b)

            #print "acb ", acb, "adb ", adb
            
            if (c_num != b_num) and ((acb > adb and acb < geo.PI) or (adb >= geo.PI and acb < geo.PI)) : # CB is edge
                self.addEdgeToPencil([c_num, b_num])

                if self.checkNewEdge(starter, [c_num, b_num]) :
                    reverse = 2
                    break
                else :
                    a_num = c_num
                    a = c
            elif (a_num != d_num) and (acb < adb and  adb < geo.PI) or (acb >= geo.PI and adb < geo.PI) : # AD is edge
                self.addEdgeToPencil([a_num, d_num])

                if self.checkNewEdge(starter, [a_num, d_num]) :
                    reverse = 2
                    break
                else :
                    b_num = d_num
                    b = d
            else :
                break
            #self.drawEdge(a, b, colorEdge, 3, 0)
            #self.drawEdge(a, b, "red", 2, 0)
            #self.drawAllPoints()
            #nb = raw_input()

        if reverse == 0 :
            self.sewTriangle([starter[1], starter[0]], 1, "blue")

    def localizePointInPencil(self, srcPoint, locatablePoint) :
        for i in range(len(self.neighbors[2][srcPoint])) :
            if self.neighbors[2][srcPoint][i] == locatablePoint :
                return i
        return -1

    def findNextStarter(self) :
        open = -1
        num_open = -1
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
                    num_open = i
                    open = num_p
                    break

        #self.drawPoint(self.points[2][open], "black", 5, 1)

        # first way
        openPoint = self.points[2][open]
        pencil_p = (open >= len(self.points[0]))
        if pencil_p == 0 :
            start = len(self.points[0])
            end = len(self.points[2])
        else :
            start = 0
            end = len(self.points[0])

        bridgeParameters = geo.lineParameters(self.points[2][curBridge[0]], self.points[2][curBridge[1]])

        minLength = 1e+10
        num_endStarter = -1;
        for i in range(start, end) :
            tempPoint = self.points[2][i]
            
            [a, b, c] = geo.lineParameters(tempPoint, openPoint)

            midlePoint = [round((tempPoint[0] + openPoint[0]) / 2), round((tempPoint[1] + openPoint[1]) / 2)]
            [pa, pb, pc] = geo.perpendicularParameters(midlePoint, [a, b, c])

            p0 = geo.intersectLines(bridgeParameters, [pa, pb, pc])
            #print p0
            #self.drawPoint(p0, "black", 5, 0)
            #nb = raw_input()

            temp = geo.getLength(p0, openPoint)
            if min(self.points[2][curBridge[0]][0], self.points[2][curBridge[1]][0]) <= p0[0] and \
               p0[0] <= max(self.points[2][curBridge[0]][0], self.points[2][curBridge[1]][0]) and \
               min(self.points[2][curBridge[0]][1], self.points[2][curBridge[1]][1]) <= p0[1] and \
               p0[1] <= max(self.points[2][curBridge[0]][1], self.points[2][curBridge[1]][1]) and \
               temp < minLength :
                minLength = temp
                num_endStarter = i

        # second way
        fix = curBridge[1 - num_open]
        fixPoint = self.points[2][fix]

        right_num = 0
        left_num = 1
        
        draw.drawEdge(self.canvas, self.points[2][curBridge[0]], self.points[2][curBridge[1]], "purple", 3)
        #self.drawPoint(self.points[2][fix], "blue", 5)

        for i in range(len(self.neighbors[2][fix])) :
            p1 = self.points[2][self.neighbors[2][fix][i]]
            p2 = self.points[2][self.neighbors[2][fix][(i + 1) % len(self.neighbors[2][fix])]]
            tempLineParameters = geo.lineParameters(p1, p2)

            p0 = geo.intersectLines(bridgeParameters, tempLineParameters)

            if geo.inSegment(fixPoint, openPoint, p0) and geo.inSegment(p1, p2, p0) :
                right_num = self.neighbors[2][fix][i]
                left_num = self.neighbors[2][fix][(i + 1) % len(self.neighbors[2][fix])]
                draw.drawEdge(self.canvas, self.points[2][fix], self.points[2][right_num], "blue", 3)
                draw.drawEdge(self.canvas, self.points[2][fix], self.points[2][left_num], "blue", 3)
                break
        
        draw.drawEdge(self.canvas, self.points[2][right_num], self.points[2][left_num], "red", 3)

        """loc1 = self.localizePointInPencil(left_num, right_num)
        loc2 = self.localizePointInPencil(left_num, fix)
        if loc1 == -1 :
            return -1
        else :
            step = geo.sign(loc1 - loc2) % len(self.neighbors[2][left_num])"""

        #print step
        #nb = raw_input()
        #print dir
        step = geo.sign(geo.exteriorProd(openPoint, fixPoint, self.points[2][right_num]))
        #print geo.exteriorProd(openPoint, fixPoint, self.points[2][left_num])
        #print "step", step
        
        while(1) :
            loc_left_in_right = self.localizePointInPencil(right_num, left_num)
            next_right = self.neighbors[2][right_num][(loc_left_in_right - step) % len(self.neighbors[2][right_num])]
            draw.drawEdge(self.canvas, self.points[2][right_num], self.points[2][next_right], "pink", 3)
                
            loc_right_in_left = self.localizePointInPencil(left_num, right_num)
            next_left = self.neighbors[2][left_num][(loc_right_in_left + step) % len(self.neighbors[2][left_num])]
            draw.drawEdge(self.canvas, self.points[2][left_num], self.points[2][next_left], "pink", 3)

            nb = raw_input()

            #print next_right_num
            tempLineParameters = geo.lineParameters(self.points[2][right_num], self.points[2][next_right])
            p0 = geo.intersectLines(bridgeParameters, tempLineParameters)
            if geo.inSegment(fixPoint, openPoint, p0) and geo.inSegment(self.points[2][right_num], self.points[2][next_right], p0) :
                draw.drawEdge(self.canvas, self.points[2][right_num], self.points[2][next_right], "red", 3)
                left_num = next_right
            else :
                tempLineParameters = geo.lineParameters(self.points[2][left_num], self.points[2][next_left])
                p0 = geo.intersectLines(bridgeParameters, tempLineParameters)
                if geo.inSegment(fixPoint, openPoint, p0) and geo.inSegment(self.points[2][left_num], self.points[2][next_left], p0) :
                    draw.drawEdge(self.canvas, self.points[2][left_num], self.points[2][next_left], "red", 3)
                    right_num = next_left
                else :
                    break
            nb = raw_input()

        draw.drawEdge(self.canvas, self.points[2][open], self.points[2][num_endStarter], "red", 3, 1)
        nb = raw_input()

        return [open, num_endStarter]

    """def checkStruct(self) :
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
            print self.points[1]"""

    def drawAllPoints(self) :
        for i in range(2) :
            for p in self.points[i] :
                draw.drawPoint(self.canvas, p, self.color[i], 3, i)

    def makeConcatenation(self, col = "black", wid = 1) :
        #print "self.points[0] =", self.points[0]
        #print "self.points[1] =", self.points[1]
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
            self.sewTriangle(starter, 0, "#00CC33")
        
        finish = time()

        self.drawStruct(col, wid)
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

        randomPoints = 0 # 0 - example, 1 - rand

        if randomPoints == 0 :
            # separeted triangle seam
            #self.points[0] = [[147, 93], [226, 45], [253, 130], [149, 197], [78, 54]]
            #self.points[1] = [[391, 124], [478, 133], [432, 181]]
            
            #example starter in report
            #self.points[0] = [[208, 180], [151, 210], [175, 126], [252, 134], [236, 227]]
            #self.points[1] = [[192, 203], [199, 147], [249, 187], [300, 158], [270, 236]]

            #model example in report
            #self.points[0] = [[139, 243], [173, 169], [237, 195], [211, 230], [275, 231], [224, 273], [178, 279], [228, 149], [286, 195], [291, 143], [333, 169], [312, 187], [355, 218], [327, 254], [293, 279], [369, 261], [394, 191], [356, 139]]
            #self.points[1] = [[187, 255], [198, 193], [248, 251], [264, 173], [315, 226], [329, 294], [262, 302], [404, 236], [358, 186], [323, 119], [411, 151], [398, 294]]
            
            #circlic seam
            #self.points[0] = [[128, 187], [158, 124], [220, 146], [220, 202], [170, 236], [173, 191]]
            #self.points[1] = [[181, 157], [226, 113], [283, 118], [273, 199], [300, 236], [237, 251]]

            # big example
            #self.points[0] = [[144, 161], [272, 247], [70, 341], [360, 205], [474, 294], [87, 153], [351, 127], [543, 136], [468, 324], [447, 38], [571, 129], [176, 182], [179, 341], [452, 212], [95, 320], [31, 236], [265, 75], [548, 21], [364, 191], [111, 78], [369, 88], [454, 105], [42, 233], [49, 125], [503, 175], [56, 93], [92, 248], [474, 92], [19, 269], [467, 283], [427, 239], [463, 15], [232, 20], [405, 141], [135, 177], [483, 37], [440, 102], [184, 160], [553, 147], [338, 330], [121, 26], [306, 152], [137, 93], [107, 124], [448, 339], [481, 77], [397, 113], [570, 285], [11, 38], [343, 113], [379, 39], [16, 146], [238, 196], [481, 125], [24, 207], [97, 81], [535, 194], [114, 46], [80, 213], [512, 268], [72, 276], [185, 21], [429, 278], [90, 212], [75, 65], [512, 7], [287, 5], [532, 284], [371, 234], [455, 311], [367, 320], [31, 151], [406, 188], [359, 252], [457, 334], [479, 185], [39, 213], [236, 298], [217, 29], [123, 275], [508, 285], [294, 67], [375, 219], [462, 118], [196, 117], [15, 179], [318, 9], [525, 290], [121, 86], [47, 207], [321, 27], [539, 206], [74, 82], [33, 43], [360, 144], [50, 101], [392, 160], [550, 126], [38, 33], [159, 110]]
            #self.points[1] = [[444, 249], [337, 338], [574, 206], [212, 173], [437, 275], [477, 118], [172, 163], [129, 232], [481, 143], [551, 65], [543, 12], [93, 88], [556, 85], [468, 236], [477, 86], [438, 159], [123, 288], [562, 14], [19, 57], [176, 256], [98, 302], [228, 343], [451, 52], [571, 211], [108, 62], [484, 243], [496, 250], [315, 342], [358, 264], [499, 39], [321, 283], [269, 228], [81, 81], [113, 245], [201, 296], [371, 250], [579, 263], [315, 105], [344, 48], [290, 57], [499, 251], [81, 238], [509, 109], [510, 246], [325, 128], [118, 302], [357, 16], [78, 217], [69, 133], [195, 343], [545, 232], [236, 231], [467, 341], [204, 329], [545, 301], [89, 166], [477, 159], [483, 49], [35, 314], [343, 327], [312, 188], [107, 244], [24, 120], [126, 165], [475, 120], [15, 72], [344, 97], [221, 44], [111, 246], [16, 182], [181, 71], [305, 267], [111, 10], [94, 187], [144, 308], [190, 105], [366, 277], [474, 231], [308, 111], [65, 17], [477, 343], [475, 81], [574, 39], [581, 36], [538, 330], [506, 167], [317, 123], [243, 213], [54, 181], [363, 113], [28, 214], [211, 337], [375, 146], [481, 279], [552, 193], [456, 270], [519, 183], [133, 20], [250, 15], [523, 95]]
            
            #self.points[0] = [[484, 230], [268, 43], [515, 237], [472, 191], [164, 285], [50, 331], [578, 200], [436, 21], [33, 276], [426, 179], [123, 13], [450, 216], [268, 218], [276, 97], [227, 331], [196, 176], [114, 338], [39, 326], [140, 76], [534, 295]]
            #self.points[1] = [[400, 201], [281, 56], [281, 185], [538, 269], [45, 253], [185, 94], [433, 37], [390, 250], [548, 23], [405, 324], [568, 298], [329, 138], [422, 257], [361, 275], [469, 344], [140, 141], [433, 264], [314, 17], [95, 114], [49, 100]]
        
            #find starter
            self.points[0] = [[103, 167], [164, 82], [273, 178], [324, 89], [406, 161]]
            self.points[1] = [[200, 143], [196, 213], [307, 132], [356, 190], [453, 97], [462, 220]]

        elif randomPoints == 1 :
            nPoints = [50, 50]
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
            #line = geo.lineParameters(p1, p2)
            #lineP = geo.perpendicular([(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2], line)
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
        #print "sciPy", sciPyTime, "my", myTime, "myTest", myTimeTest
        #print "my", myTime, "myTest", myTimeTest

    def experimentTime(self) :
        fout = open('res.txt', 'w')
        nPoints = [10 * x for x in range(1, 71)]
        for n in nPoints :
            sciPyTime = 0
            myTime = 0
            numIter = 10
            for it in range(numIter) :
                while True :
                    self.startTriDone = False
                    self.checkTriDone = False
                    self.firstTri = True
                    self.points = [[], [], []]
                    self.faces = [[], [], []]
                    self.bridges = []

                    self.canvas.delete("all")
                    self.experimentMode = True

                    for i in range(2) :
                        x = np.random.randint(0, self.width - 20, (n, 1)) + 10
                        y = np.random.randint(0, self.height - 60, (n, 1)) + 5
                        self.points[i] = np.concatenate((x, y), axis = 1)

                        self.secondTriangle()

                    self.points[0] = self.points[0].tolist()
                    self.points[1] = self.points[1].tolist()


                    try :
                        [sciPy, my] = self.makeTriangle()
                        sciPyTime += sciPy
                        myTime += my
                        break
                    except :
                        print "exception"
                        pass
            print n
            fout.write(str(n) + ' ' + str(sciPyTime / numIter) + ' ' + str(myTime / numIter) + '\n')