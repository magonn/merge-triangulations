from Tkinter import *
import draw

from random import *
from time import *

from scipy.spatial import Delaunay
import numpy as np
import Queue

import geo

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
            
            for i in xrange(2) :
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
        for j in xrange(len(self.tree[num])) :
            p1 = self.points[num][self.tree[num][j][0]]
            p2 = self.points[num][self.tree[num][j][1]]
            
            draw.drawEdge(self.canvas, p1, p2, col, wid)

    def makeMST(self) :
        edgesMST = [[], []]
        self.tree = [[], []]
        
        for i in xrange(2) :
            for f in self.faces[i] :
                p = [f[0], f[1], f[2], f[0]]
            
                for j in xrange(3) :
                    lenEdge = geo.getLength(self.points[i][p[j]], self.points[i][p[j + 1]])

                    if edgesMST[i].count((lenEdge, p[j], p[j + 1])) == 0 :
                        edgesMST[i].append((lenEdge, p[j], p[j + 1]))

            edgesMST[i].sort()
            
            treeId = []
            for j in xrange(len(self.points[i])) :
                treeId.append(j)

            for e in edgesMST[i] :
                a = e[1]
                b = e[2]

                if treeId[a] != treeId[b] :
                    self.tree[i].append((a, b))
                    oldId = treeId[b]
                    newId = treeId[a]
                    for j in xrange(len(treeId)) :
                        if treeId[j] == oldId :
                            treeId[j] = newId

            #self.drawMST(i, "black", 2)

        self.mst = []
        n1 = len(self.points[0])
        
        for i in xrange(2) :
            for j in xrange(len(self.tree[i])) :
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
        for i in xrange(2) :
            for j in xrange(len(self.points[i])) :
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
        for j in xrange(len(self.points[1 - fixPoint])) :
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
        for i in xrange(len(self.points[2])) :
            p1 = self.points[2][i]
            for j in self.neighbors[2][i] :
                draw.drawEdge(self.canvas, p1, self.points[2][j], col, wid, 0)

    def createStruct(self) :
        self.neighbors = [[], [], []]

        for i in xrange(2) :
            self.nPoints[i] = len(self.points[i])
            temp = [[] for x in xrange(self.nPoints[i])]
            flag = [False for x in xrange(len(self.faces[i]))]
            
            queue = Queue.Queue()
            queue.put(0)

            while not queue.empty() :
                cur = queue.get()
                if flag[cur] == True :
                    continue

                flag[cur] = True
                
                for j in xrange(3) :
                    if self.neighFaces[i][cur][j] != -1 and flag[self.neighFaces[i][cur][j]] == False:
                        queue.put(self.neighFaces[i][cur][j])

                x = self.faces[i][cur]
                tempFace = [x[0], x[1], x[2], x[0], x[1]]
                
                for j in xrange(1, 4) :
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
        self.neighbors[2] = [[] for i in xrange(self.nPoints[2])]

        for i in xrange(2) :
            for j in xrange(self.nPoints[i]) :
                for k in xrange(len(self.neighbors[i][j])) :
                    self.neighbors[2][i * self.nPoints[0] + j].append(self.neighbors[i][j][k] + i * self.nPoints[0])

    def addPointToPencil(self, where, pNum) :
        if self.neighbors[2][where].count(pNum) > 0 :
            return

        numPoints = len(self.neighbors[2][where])
        resAngle = []
        for j in xrange(numPoints) :
            res = geo.countAngle(self.points[2][self.neighbors[2][where][j]], self.points[2][where], self.points[2][pNum])
            resAngle.append(res)
        
        if len(resAngle) == 0 :
            aMin = 0
        else :
            aMin = resAngle.index(min(resAngle))
        
        self.neighbors[2][where].insert(aMin, pNum)
        
    def addEdgeToPencil(self, edge) :
        for i in xrange(2) :
            self.addPointToPencil(edge[i], edge[1 - i])

    def checkBridge(self, edge) :
        if self.mst.count(edge) > 0 or self.mst.count([edge[1], edge[0]]) > 0 :
            self.bridges.append(edge)
        
    def deleteWrongEdges(self, starter, reverse = 0) :
        [aNum, bNum] = starter
        a = self.points[2][aNum]
        b = self.points[2][bNum]
        while(1) :
            idx = self.neighbors[2][aNum].index(bNum)
            numPoints = len(self.neighbors[2][aNum])
            c1Num = self.neighbors[2][aNum][(idx - 1) % numPoints]
            c1 = self.points[2][c1Num]
            c2Num = self.neighbors[2][aNum][(idx - 2) % numPoints]
            c2 = self.points[2][c2Num]
            
            if c1Num == bNum or c2Num == bNum :
                break
                
            abc1 = geo.countAngle(c1, b, a)
            ac2c1 = geo.countAngle(a, c2, c1)
            
            if abc1 + ac2c1 <= geo.PI or abc1 > geo.PI or ac2c1 > geo.PI:
                break
            else :
                self.checkBridge([aNum, c1Num])
                self.neighbors[2][aNum].remove(c1Num)
                self.neighbors[2][c1Num].remove(aNum)
                
        while(1) :
            idx = self.neighbors[2][aNum].index(bNum)
            numPoints = len(self.neighbors[2][aNum])
            c1Num = self.neighbors[2][aNum][(idx + 1) % numPoints]
            c1 = self.points[2][c1Num]
            c2Num = self.neighbors[2][aNum][(idx + 2) % numPoints]
            c2 = self.points[2][c2Num]

            if c1Num == bNum or c2Num == bNum :
                break
                
            abc1 = geo.countAngle(a, b, c1)
            ac2c1 = geo.countAngle(c1, c2, a)
            
            if abc1 + ac2c1 <= geo.PI or abc1 > geo.PI or ac2c1 > geo.PI :
                break
            else :
                self.checkBridge([aNum, c1Num])
                self.neighbors[2][aNum].remove(c1Num)
                self.neighbors[2][c1Num].remove(aNum)
                
        if reverse == 0 :
            self.deleteWrongEdges([starter[1], starter[0]], 1)

    def checkNewEdge(self, e1, e2) :
        return (e1[0] == e2[0] and e1[1] == e2[1]) or (e1[0] == e2[1] and e1[1] == e2[0])

    def sewTriangle(self, starter, reverse = 0, colorEdge = "purple") :
        [bNum, aNum] = starter
        a = self.points[2][aNum]
        b = self.points[2][bNum]
        #self.drawEdge(a, b, colorEdge, 2)
        #self.drawEdge(a, b, "black", 3)
        #self.drawAllPoints()
        #nb = raw_input()
        
        while(1) :
            self.deleteWrongEdges([aNum, bNum])
            
            idx1 = self.neighbors[2][aNum].index(bNum)
            cNum = self.neighbors[2][aNum][(idx1 - 1) % len(self.neighbors[2][aNum])]
            c = self.points[2][cNum]

            idx2 = self.neighbors[2][bNum].index(aNum)
            dNum = self.neighbors[2][bNum][(idx2 + 1) % len(self.neighbors[2][bNum])]
            d = self.points[2][dNum]
            
            #print "EDGE", aNum, bNum, cNum, dNum
            #self.drawPoint(c, "yellow", 3, 0)
            #self.drawPoint(d, "white", 3, 0)
            #nb = raw_input()

            acb = geo.countAngle(a, c, b)
            adb = geo.countAngle(a, d, b)

            #print "acb ", acb, "adb ", adb
            
            if (cNum != bNum) and ((acb > adb and acb < geo.PI) or (adb >= geo.PI and acb < geo.PI)) : # CB is edge
                self.addEdgeToPencil([cNum, bNum])

                if self.checkNewEdge(starter, [cNum, bNum]) :
                    reverse = 2
                    break
                else :
                    aNum = cNum
                    a = c
            elif (aNum != dNum) and (acb < adb and  adb < geo.PI) or (acb >= geo.PI and adb < geo.PI) : # AD is edge
                self.addEdgeToPencil([aNum, dNum])

                if self.checkNewEdge(starter, [aNum, dNum]) :
                    reverse = 2
                    break
                else :
                    bNum = dNum
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
        for i in xrange(len(self.neighbors[2][srcPoint])) :
            if self.neighbors[2][srcPoint][i] == locatablePoint :
                return i
        return -1

    def findNextStarter(self) :
        [fix, open] = self.bridges.pop(0)

        ### define open and fix points
        openPoint = self.points[2][open]
        fixPoint = self.points[2][fix]
        bridgeParameters = geo.lineParameters(openPoint, fixPoint)

        ### first way
        """if open < len(self.points[0]) :
            start = len(self.points[0])
            end = len(self.points[2])
        else :
            start = 0
            end = len(self.points[0])

        minRadius = geo.getLength(openPoint, fixPoint) / 2
        endStarterNum = -1;
        for i in xrange(start, end) :
            tempPoint = self.points[2][i]
            
            radius = geo.radiusByLineAndPoint(bridgeParameters, openPoint, fixPoint, tempPoint)
            
            if radius != -1 and radius < minRadius :
                minRadius = radius
                endStarterNum = i"""

        ### second way
        rightNum = 0
        leftNum = 1
        
        draw.drawEdge(self.canvas, openPoint, fixPoint, "purple", 3) # draw bridge
        
        ### find left and right neighbors of bridge
        for i in xrange(len(self.neighbors[2][fix])) :
            p1Num = self.neighbors[2][fix][i]
            p1 = self.points[2][p1Num]
            p2Num = self.neighbors[2][fix][(i + 1) % len(self.neighbors[2][fix])]
            p2 = self.points[2][p2Num]
            
            lineParameters = geo.lineParameters(p1, p2)

            p0 = geo.intersectLines(bridgeParameters, lineParameters)

            if geo.inSegment(fixPoint, openPoint, p0) and geo.inSegment(p1, p2, p0) and geo.countAngle(p2, fixPoint, p1) < geo.pi:
                rightNum = p1Num
                leftNum = p2Num
                rightPoint = p1
                leftPoint = p2
                #draw.drawEdge(self.canvas, self.points[2][fix], self.points[2][rightNum], "blue", 3)
                #draw.drawEdge(self.canvas, self.points[2][fix], self.points[2][leftNum], "blue", 3)
                break
        
        #draw.drawEdge(self.canvas, self.points[2][rightNum], self.points[2][leftNum], "red", 3) # draw intersection line
        #nb = raw_input()

        step = geo.sign(geo.exteriorProd(openPoint, fixPoint, self.points[2][rightNum]))
        
        minRadiusTemp = geo.getLength(openPoint, fixPoint) / 2
        endStarterNumTemp = -1

        ### update minRadius with left and right neighbours
        radius = geo.radiusByLineAndPoint(bridgeParameters, openPoint, fixPoint, rightPoint)
        if radius != -1 and radius < minRadiusTemp :
            minRadiusTemp = radius
            endStarterNumTemp = rightNum

        radius = geo.radiusByLineAndPoint(bridgeParameters, openPoint, fixPoint, leftPoint)
        if radius != -1 and radius < minRadiusTemp :
            minRadiusTemp = radius
            endStarterNumTemp = leftNum

        checkPoints = []
        
        while(1) :
            ### find next intersection
            locLeftInRight = self.localizePointInPencil(rightNum, leftNum)
            nextRightNum = self.neighbors[2][rightNum][(locLeftInRight - step) % len(self.neighbors[2][rightNum])] 
            
            locRightInLeft = self.localizePointInPencil(leftNum, rightNum)
            nextLeftNum = self.neighbors[2][leftNum][(locRightInLeft + step) % len(self.neighbors[2][leftNum])]
            
            #draw.drawEdge(self.canvas, self.points[2][rightNum], self.points[2][nextRightNum], "pink", 3)
            #draw.drawEdge(self.canvas, self.points[2][leftNum], self.points[2][nextLeftNum], "pink", 3)
            
            if nextRightNum != nextLeftNum :
                checkPoints.extend(self.neighbors[2][leftNum])
                checkPoints.extend(self.neighbors[2][rightNum])
                break

            nextNum = nextRightNum
            nextPoint = self.points[2][nextNum]

            #nb = raw_input()

            lineParameters = geo.lineParameters(rightPoint, nextPoint)
            p0 = geo.intersectLines(bridgeParameters, lineParameters)

            if geo.inSegment(fixPoint, openPoint, p0) and geo.inSegment(rightPoint, nextPoint, p0) and \
                p0 != rightPoint and p0 != nextPoint :
                #draw.drawEdge(self.canvas, rightPoint, nextPoint, "red", 3) # draw intersection line
                leftNum = nextNum
                leftPoint = nextPoint

                radius = geo.radiusByLineAndPoint(bridgeParameters, openPoint, fixPoint, leftPoint)
                if radius != -1 and radius < minRadiusTemp :
                    minRadiusTemp = radius
                    endStarterNumTemp = leftNum
            else :
                lineParameters = geo.lineParameters(leftPoint, nextPoint)
                p0 = geo.intersectLines(bridgeParameters, lineParameters)
                
                if geo.inSegment(fixPoint, openPoint, p0) and geo.inSegment(leftPoint, nextPoint, p0) and \
                    p0 != leftPoint and p0 != nextPoint :
                    #draw.drawEdge(self.canvas, leftPoint, nextPoint, "red", 3) # draw intersection line
                    rightNum = nextNum
                    rightPoint = nextPoint

                    radius = geo.radiusByLineAndPoint(bridgeParameters, openPoint, fixPoint, rightPoint)
                    if radius != -1 and radius < minRadiusTemp :
                        minRadiusTemp = radius
                        endStarterNumTemp = rightNum
                else :
                    checkPoints.extend(self.neighbors[2][leftNum])
                    checkPoints.extend(self.neighbors[2][rightNum])
                    checkPoints.extend(self.neighbors[2][nextNum])
                    break

            #nb = raw_input()
        
        for i in xrange(len(checkPoints)) :
            p = self.points[2][checkPoints[i]]
            radius = geo.radiusByLineAndPoint(bridgeParameters, openPoint, fixPoint, p)
            if radius != -1 and radius < minRadiusTemp :
                minRadiusTemp = radius
                endStarterNumTemp = checkPoints[i]

        #draw.drawEdge(self.canvas, openPoint, self.points[2][endStarterNum], "yellow", 3) # next draw starter
        draw.drawEdge(self.canvas, openPoint, self.points[2][endStarterNumTemp], "black", 3) # next draw starter
        #nb = raw_input()
        
        """if endStarterNum != endStarterNumTemp :
            print "FUUUUUUUUUUU"
            print minRadius, minRadiusTemp
            print endStarterNum, endStarterNumTemp
            print "self.points[0] =", self.points[0]
            print "self.points[1] =", self.points[1]"""      

        return [open, endStarterNumTemp]

    """def checkStruct(self) :
        flag = False
        for i in xrange(len(self.points[2])) :
            for j in xrange(len(self.neighbors[2][i])) :
                pNum = self.neighbors[2][i][j]
                if self.neighbors[2][pNum].count(i) != 1 :
                    flag = True
                    print "bad struct"
                    print i, self.neighbors[2][i]
                    print pNum, self.neighbors[2][pNum]
                    break

        if flag :
            print "There are errors"
            print self.points[0]
            print self.points[1]"""

    def drawAllPoints(self) :
        for i in xrange(2) :
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
            #print "find starter ...."
            starter = self.findNextStarter()
            #print "starter"
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

        randomPoints = 1 # 0 - example, 1 - rand

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

            self.points[0] = [[457, 272], [90, 144], [321, 68], [383, 303], [75, 245], [315, 47], [109, 65], [197, 298], [376, 241], [441, 37], [298, 331], [436, 233], [122, 206], [224, 339], [553, 198], [38, 65], [470, 245], [27, 9], [183, 233], [519, 241], [286, 203], [283, 132], [443, 196], [458, 64], [233, 71], [254, 269], [251, 72], [409, 150], [493, 33], [313, 21], [24, 283], [105, 62], [537, 327], [318, 292], [217, 21], [40, 201], [425, 257], [11, 292], [89, 183], [341, 148], [149, 66], [442, 23], [565, 195], [132, 115], [22, 212], [384, 67], [101, 187], [576, 249], [472, 234], [128, 340]]
            self.points[1] = [[168, 157], [208, 337], [472, 234], [24, 39], [202, 230], [552, 173], [337, 44], [256, 47], [114, 124], [33, 37], [114, 103], [538, 93], [265, 7], [329, 333], [147, 108], [491, 334], [260, 240], [394, 182], [269, 186], [414, 241], [578, 282], [149, 292], [448, 57], [219, 80], [202, 316], [31, 206], [63, 268], [186, 162], [329, 85], [476, 199], [141, 257], [113, 174], [273, 30], [398, 5], [423, 340], [226, 49], [95, 237], [406, 12], [57, 63], [58, 27], [380, 7], [208, 21], [469, 102], [466, 57], [97, 225], [116, 194], [180, 117], [216, 50], [115, 305], [86, 127]]
        elif randomPoints == 1 :
            nPoints = [50, 50]
            for i in xrange(2) :
                x = np.random.randint(0, self.width - 20, (nPoints[i], 1)) + 10
                y = np.random.randint(0, self.height - 60, (nPoints[i], 1)) + 5
                self.points[i] = np.concatenate((x, y), axis = 1)

                self.secondTriangle()

            self.points[0] = self.points[0].tolist()
            self.points[1] = self.points[1].tolist()
            print "self.points[0] =", self.points[0]
            print "self.points[1] =", self.points[1]
        else :
            p1 = [400, 200]
            p2 = [300, 250]
            draw.drawEdge(self.canvas, [200, 200], [400, 200], "black", 1, 0)
            draw.drawEdge(self.canvas, [300, 250], [400, 200], "black", 1, 1)
            draw.drawPoint(self.canvas, p1, "black", 3, 0)
            draw.drawPoint(self.canvas, p2, "black", 3, 1)
            x = 320
            y = 2 * x - 475
            #line = geo.lineParameters(p1, p2)
            #lineP = geo.perpendicular([(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2], line)
            draw.drawEdge(self.canvas, [x, y], [(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2], "black", 1, 1)
            
            draw.drawPoint(self.canvas, [(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2], "black", 3, 1)
            draw.drawPoint(self.canvas, [337, 200], "black", 3, 1)
            draw.drawPoint(self.canvas, [250, 150], "black", 3, 1)
            draw.drawPoint(self.canvas, [230, 220], "black", 3, 1)
            draw.drawPoint(self.canvas, [380, 100], "black", 3, 1)
            draw.drawCircle(self.canvas, [337, 200], "black", 400 - 337, 0)
            
            return

        #self.drawAllPoints()
        [sciPyTime, myTime] = self.makeTriangle()
        #print "sciPy", sciPyTime, "my", myTime, "myTest", myTimeTest
        #print "my", myTime, "myTest", myTimeTest

    def experimentTime(self) :
        fout = open('res.txt', 'w')
        nPoints = [10 * x for x in xrange(1, 71)]
        for n in nPoints :
            sciPyTime = 0
            myTime = 0
            numIter = 10
            for it in xrange(numIter) :
                while True :
                    self.startTriDone = False
                    self.checkTriDone = False
                    self.firstTri = True
                    self.points = [[], [], []]
                    self.faces = [[], [], []]
                    self.bridges = []

                    self.canvas.delete("all")
                    self.experimentMode = True

                    for i in xrange(2) :
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