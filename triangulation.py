from form import BaseForm

from scipy.spatial import Delaunay
import numpy as np
import Queue
import sets
from time import *

import geo
import draw

class ConstructTriangulation(BaseForm):
    
    preprocessingDone = False

    # in experement mode do not remove last point (click on button)
    experimentMode = False

    # data for creation structure
    faces = [[], [], []]
    neighFaces = [[], [], []]
    answer = []

    # data for detecting next starter
    neighbors = [[], [], []]
    fictiveEdges = []
    
    scipyTime = -1

    def Reset(self):
        self.canvas.delete("all")

        self.firstPointsSet = True
        self.preprocessingDone = False
        self.experimentMode = True
        
        self.points = [[], [], []]

        self.faces = [[], [], []]
        self.neighFaces = [[], [], []]
        self.answer = []

        self.neighbors = [[], [], []]
        self.fictiveEdges = []
        
        self.scipyTime = -1

    def Preprocessing(self):
    	if not self.preprocessingDone:
            for i in xrange(2):
                if not self.experimentMode:
                    self.points[i].pop(len(self.points[i]) - 1)
                
                tri = Delaunay(self.points[i])
                self.faces[i] = tri.vertices
                self.neighFaces[i] = tri.neighbors
                #draw.Triangle(self.canvas, self.points[i], self.faces[i], "black", 1, i)
                
            self.points[2].extend(self.points[0])
            self.points[2].extend(self.points[1])

            start = time()
            tri = Delaunay(self.points[2])
            self.scipyTime = time() - start
            
            self.faces[2] = tri.vertices
            self.neighFaces[2] = tri.neighbors
            
            #draw.Triangle(self.canvas, self.points[2], self.faces[2], "green", 3, 2)
            
        self.CreateStruct()
        self.preprocessingDone = True
            
    def DrawStruct(self, col = "black", wid = 1):
        for i in xrange(len(self.points[2])):
            p1 = self.points[2][i]
            for j in self.neighbors[2][i]:
                draw.Edge(self.canvas, p1, self.points[2][j], col, wid, 0)

    def CreateStruct(self):
        if not self.preprocessingDone:
            for i in xrange(3):
                temp = [[] for x in xrange(len(self.points[i]))]
                flag = [False for x in xrange(len(self.faces[i]))]
            
                queue = Queue.Queue()
                queue.put(0)

                while not queue.empty():
                    cur = queue.get()
                    if flag[cur] == True:
                        continue

                    flag[cur] = True
                    
                    for j in xrange(3):
                        if self.neighFaces[i][cur][j] != -1 and flag[self.neighFaces[i][cur][j]] == False:
                            queue.put(self.neighFaces[i][cur][j])

                    x = self.faces[i][cur]
                    tempFace = [x[0], x[1], x[2], x[0], x[1]]
                    
                    for j in xrange(1, 4):
                        n = len(temp[tempFace[j]])
                        if n == 0:
                            angle = geo.ExteriorProd(self.points[i][tempFace[j - 1]], self.points[i][tempFace[j]], self.points[i][tempFace[j + 1]])
                            if angle < 0:
                                temp[tempFace[j]] = [tempFace[j - 1], tempFace[j + 1]]
                            else:
                                temp[tempFace[j]] = [tempFace[j + 1], tempFace[j - 1]]
                        else:
                            c1 = temp[tempFace[j]].count(tempFace[j - 1])
                            c2 = temp[tempFace[j]].count(tempFace[j + 1])
                            
                            if temp[tempFace[j]][0] == tempFace[j - 1] and c2 == 0:
                                temp[tempFace[j]].insert(0, tempFace[j + 1])
                                c2 = 1

                            elif temp[tempFace[j]][0] == tempFace[j + 1] and c1 == 0:
                                temp[tempFace[j]].insert(0, tempFace[j - 1])
                                c1 = 1

                            elif temp[tempFace[j]][n - 1] == tempFace[j - 1] and c2 == 0:
                                temp[tempFace[j]].insert(n, tempFace[j + 1])
                                c2 = 1

                            elif temp[tempFace[j]][n - 1] == tempFace[j + 1] and c1 == 0:
                                temp[tempFace[j]].insert(n, tempFace[j - 1])
                                c1 = 1

                            if c1 == 0 or c2 == 0:
                                flag[cur] = False
                                queue.put(cur)

                self.neighbors[i] = temp

            self.answer = self.neighbors[2]

        # save only initial struct
        self.neighbors[2] = [[] for i in xrange(len(self.points[2]))]

        for i in xrange(2):
            for j in xrange(len(self.points[i])):
                for k in xrange(len(self.neighbors[i][j])):
                    self.neighbors[2][i * len(self.points[0]) + j].append(self.neighbors[i][j][k] + i * len(self.points[0]))

    def GetFirstStarter(self):
        twoStarter = [0, 0]
        for i in xrange(2):
            for j in xrange(len(self.points[i])):
                p = self.points[i][j]
                tempStarter = self.points[i][twoStarter[i]]
                if (p[0] < tempStarter[0]) or (p[0] == tempStarter[0] and p[1] > tempStarter[1]):
                    twoStarter[i] = j
        
        pTwoStarter = [self.points[0][twoStarter[0]], self.points[1][twoStarter[1]]]
        fixPoint = (pTwoStarter[0][0] < pTwoStarter[1][0]) or (pTwoStarter[0][0] == pTwoStarter[1][0] and (pTwoStarter[0][1] > pTwoStarter[1][1]))
        
        pEndFixPoint = pTwoStarter[1 - fixPoint]
        endFixPoint = twoStarter[1 - fixPoint]
        
        minLength = 1e+10
        for j in xrange(len(self.points[1 - fixPoint])):
            p = self.points[1 - fixPoint][j]
            if p[0] == pTwoStarter[fixPoint][0]:
                continue
            
            temp = ((p[0] - pTwoStarter[fixPoint][0]) ** 2 + (p[1] - pTwoStarter[fixPoint][1]) ** 2) / (2 * (pTwoStarter[fixPoint][0] - p[0]))
            if 0 < temp and temp < minLength:
                minLength = temp
                endFixPoint = j
                pEndFixPoint = p
                
        pRes = [pTwoStarter[fixPoint], pEndFixPoint]
        if fixPoint == 1:
            res = [endFixPoint, twoStarter[1]]
        else:
            res = [twoStarter[0], endFixPoint]
        res = [res[0], res[1] + len(self.points[0])]
        return res

    def AddPointToPencil(self, where, pNum):
        if self.neighbors[2][where].count(pNum) > 0:
            return False

        numPoints = len(self.neighbors[2][where])
        resAngle = []
        for j in xrange(numPoints):
            res = geo.CountAngle(self.points[2][self.neighbors[2][where][j]], self.points[2][where], self.points[2][pNum])
            resAngle.append(res)
        
        if len(resAngle) == 0:
            aMin = 0
        else:
            aMin = resAngle.index(min(resAngle))
        
        self.neighbors[2][where].insert(aMin, pNum)
        return True
        
    def AddEdgeToPencil(self, edge):
        return self.AddPointToPencil(edge[0], edge[1]) and self.AddPointToPencil(edge[1], edge[0])

    def IsFictiveEdge(self, edge, breaker):
        [aNum, c1Num] = edge
        
        num1 = int(aNum >= len(self.points[0]))
        num2 = int(c1Num >= len(self.points[0]))
        
        if num1 != num2:
            print "Noooooooooooooo"
            return
        num = num1
        
        aNum -= num * len(self.points[0])
        c1Num -= num * len(self.points[0])

        idx = self.neighbors[num][aNum].index(c1Num)
        numPoints = len(self.neighbors[num][aNum])
        c2Num = self.neighbors[num][aNum][(idx - 1) % numPoints]
        dNum = self.neighbors[num][aNum][(idx + 1) % numPoints]

        adc1 = geo.CountAngle(self.points[num][c1Num], self.points[num][dNum], self.points[num][aNum])
        ac2c1 = geo.CountAngle(self.points[num][aNum], self.points[num][c2Num], self.points[num][c1Num])

        if (adc1 <= geo.PI / 2  or adc1 > geo.PI) and (ac2c1 <= geo.PI / 2 or ac2c1 > geo.PI):
            self.fictiveEdges.append([edge[0], edge[1], breaker])
            
    def DeleteFictiveEdges(self, starter, reverse = 0):
        [aNum, bNum] = starter
        a = self.points[2][aNum]
        b = self.points[2][bNum]
        
        while(1):
            idx = self.neighbors[2][aNum].index(bNum)
            numPoints = len(self.neighbors[2][aNum])
            c1Num = self.neighbors[2][aNum][(idx - 1) % numPoints]
            c1 = self.points[2][c1Num]
            c2Num = self.neighbors[2][aNum][(idx - 2) % numPoints]
            c2 = self.points[2][c2Num]
            
            if c1Num == bNum or c2Num == bNum:
                break
                
            abc1 = geo.CountAngle(c1, b, a)
            ac2c1 = geo.CountAngle(a, c2, c1)
            
            if abc1 + ac2c1 <= geo.PI or abc1 > geo.PI or ac2c1 > geo.PI:
                break
            else:
                self.IsFictiveEdge([aNum, c1Num], bNum)

                self.neighbors[2][aNum].remove(c1Num)
                self.neighbors[2][c1Num].remove(aNum)
                
        while(1):
            idx = self.neighbors[2][aNum].index(bNum)
            numPoints = len(self.neighbors[2][aNum])
            c1Num = self.neighbors[2][aNum][(idx + 1) % numPoints]
            c1 = self.points[2][c1Num]
            c2Num = self.neighbors[2][aNum][(idx + 2) % numPoints]
            c2 = self.points[2][c2Num]

            if c1Num == bNum or c2Num == bNum:
                break
                
            abc1 = geo.CountAngle(a, b, c1)
            ac2c1 = geo.CountAngle(c1, c2, a)
            
            if abc1 + ac2c1 <= geo.PI or abc1 > geo.PI or ac2c1 > geo.PI:
                break
            else:
                self.IsFictiveEdge([aNum, c1Num], bNum)
                    
                self.neighbors[2][aNum].remove(c1Num)
                self.neighbors[2][c1Num].remove(aNum)
                
        if reverse == 0:
            self.DeleteFictiveEdges([starter[1], starter[0]], 1)

    def IsEdgeExists(self, e1, e2):
        return (e1[0] == e2[0] and e1[1] == e2[1]) or (e1[0] == e2[1] and e1[1] == e2[0])

    def SewTriangles(self, starter, reverse = 0, colorEdge = "purple"):
        [aNum, bNum] = starter
        a = self.points[2][aNum]
        b = self.points[2][bNum]
            
        #draw.Edge(self.canvas, self.points[2][aNum], self.points[2][bNum], colorEdge, 2)
        #raw_input()

        while(1):
            self.DeleteFictiveEdges([aNum, bNum])
            
            idx1 = self.neighbors[2][aNum].index(bNum)
            cNum = self.neighbors[2][aNum][(idx1 - 1) % len(self.neighbors[2][aNum])]
            c = self.points[2][cNum]

            idx2 = self.neighbors[2][bNum].index(aNum)
            dNum = self.neighbors[2][bNum][(idx2 + 1) % len(self.neighbors[2][bNum])]
            d = self.points[2][dNum]
            
            acb = geo.CountAngle(a, c, b)
            adb = geo.CountAngle(a, d, b)

            if (cNum != bNum) and ((acb > adb and acb < geo.PI) or (adb >= geo.PI and acb < geo.PI)): # CB is edge
                #draw.Edge(self.canvas, self.points[2][cNum], self.points[2][bNum], colorEdge, 2)
                #raw_input()

                self.AddEdgeToPencil([cNum, bNum])

                if self.IsEdgeExists(starter, [cNum, bNum]):
                    reverse = 2
                    break
                else:
                    aNum = cNum
                    a = c
            elif (aNum != dNum) and (acb < adb and  adb < geo.PI) or (acb >= geo.PI and adb < geo.PI): # AD is edge
                #draw.Edge(self.canvas, self.points[2][aNum], self.points[2][dNum], colorEdge, 2)
                #raw_input()

                self.AddEdgeToPencil([aNum, dNum])

                if self.IsEdgeExists(starter, [aNum, dNum]):
                    reverse = 2
                    break
                else:
                    bNum = dNum
                    b = d
            else:
                break
            
        if reverse == 0:
            self.SewTriangles([starter[1], starter[0]], 1, "blue")

    def GetNextStarter(self):
        [fix, open, breaker] = self.fictiveEdges.pop(0)
        
        center = [(self.points[2][open][0] + self.points[2][fix][0]) / 2, (self.points[2][open][1] + self.points[2][fix][1]) / 2]
        
        openPoint = self.points[2][open]
        fixPoint = self.points[2][fix]

        radiusBridge = geo.GetLength(center, openPoint)
        bridgeParameters = geo.LineParameters(openPoint, center)
        
        #draw.Edge(self.canvas, self.points[2][fix], self.points[2][open], "red", 3)
        #draw.Point(self.canvas, self.points[2][breaker], "green", 5)
        #draw.Circle(self.canvas, center, "blue", wid = radiusBridge)
        #draw.AllPoints(self.canvas, self.points)
        #raw_input()

        checkPoints = Queue.Queue()
        checkPoints.put(breaker)

        visited = set()

        minRadius = -1
        endStarterNum = -1

        while not checkPoints.empty():
            nP = checkPoints.get()
            p = self.points[2][nP]

            if nP not in visited:
                visited.add(nP)
                radius = geo.RadiusByLineAndPoint(bridgeParameters, radiusBridge, center, openPoint, p)
                if radius != -1 and (minRadius == -1 or radius < minRadius):
                    minRadius = radius
                    endStarterNum = nP

                for test in self.neighbors[2][nP]:
                    if geo.GetLength(center, self.points[2][test]) < radiusBridge:
                        checkPoints.put(test)

        #draw.Edge(self.canvas, self.points[2][open], self.points[2][endStarterNum], "purple", 3)
        return [open, endStarterNum]

    def MergeTriangles(self):
        #draw.AllPoints(self.canvas, self.points)
        start = time()   

        starter = self.GetFirstStarter()
        self.AddEdgeToPencil(starter)

        self.SewTriangles(starter)
        
        all_starters = 0
        used_starters = 0
        while(len(self.fictiveEdges)):
            starter = self.GetNextStarter()
            all_starters += 1

            if self.AddEdgeToPencil(starter):
                used_starters += 1
                self.SewTriangles(starter, 0, "#00CC33")
            
        fictiveTime = time() - start

        #draw.Triangle(self.canvas, self.points[2], self.faces[2], "green", 3, 2)
        self.DrawStruct("black", 1)
        draw.AllPoints(self.canvas, self.points)

        self.CheckTriangle()
        print "ok; n =", len(self.points[2]), used_starters, "/", all_starters, "time:", fictiveTime
        return fictiveTime
    
    def CheckTriangle(self):
        for i in xrange(len(self.neighbors[2])):
            fictiveRes = self.neighbors[2][i]
            rightRes = self.answer[i]

            n = len(fictiveRes)
            if (n != len(rightRes)):
                draw.Circle(self.canvas, self.points[2][i], "green4", 5, 0)
                print "### cur", i, "res", fictiveRes, "ans", rightRes
                raise Exception("wrong output len")
                
            idx = rightRes.index(fictiveRes[0])
            for j in xrange(n):
                if fictiveRes[j] != rightRes[(idx + j) % n]:
                    draw.Circle(self.canvas, self.points[2][i], "VioletRed1", 5, 0)
                    print "res", fictiveRes, "ans", rightRes
                    raise Exception("wrong output struct")

    def Run(self):
        self.Preprocessing()
        fictiveTime = self.MergeTriangles()
        return fictiveTime

    def IncreaseDistance(self, data):
        for i in xrange(len(data) - 1):
            for j in xrange(i + 1, len(data)):
                cur = geo.GetLength(data[i], data[j])
                if cur <= geo.MINDIST:
                    data[i] = []
                    break
                    
        data = [x for x in data if len(x)]
        return data

    def GenerateData(self, n):
        x = np.random.randint(0, self.width - 20, (n, 1)) + 10
        y = np.random.randint(0, self.height - 60, (n, 1)) + 5
        data = np.concatenate((x, y), axis = 1)
        data = self.IncreaseDistance(data.tolist())
        
        delim = int(len(data) / 2)
        self.points[0] = data[:delim]
        self.points[1] = data[delim:]
        
    def Experiment(self):
        self.Reset()
        
        randomPoints = 0 # 0 - example, 1 - rand

        if randomPoints == 0:
            # model example in report
            #self.points[0] = [[139, 243], [173, 169], [237, 195], [211, 230], [275, 231], [224, 273], [178, 279], [228, 149], [286, 195], [291, 143], [333, 169], [312, 187], [355, 218], [327, 254], [293, 279], [369, 261], [394, 191], [356, 139]]
            #self.points[1] = [[187, 255], [198, 193], [248, 251], [264, 173], [315, 226], [329, 294], [262, 302], [404, 236], [358, 186], [323, 119], [411, 151], [398, 294]]
            
            # circlic seam
            #self.points[0] = [[128, 187], [158, 124], [220, 146], [220, 202], [170, 236], [173, 191]]
            #self.points[1] = [[181, 157], [226, 113], [283, 118], [273, 199], [300, 236], [237, 251]]

            # Error
            #self.points[0] = [[200, 142], [384, 182], [467, 277], [15, 160], [234, 147]]
            #self.points[1] = [[201, 270], [509, 277], [315, 241], [262, 79], [54, 277]]
            self.points[0] = [[204, 191], [397, 243], [348, 130], [145, 80], [501, 24]]
            self.points[1] = [[495, 24], [194, 89], [443, 233], [211, 239], [270, 60]]

        else:
            n = 100
            self.GenerateData(n)

        try:
            self.Run()
        except Exception as error:
            print "!!! ERROR !!!"
            print "self.points[0] =", self.points[0]
            print "self.points[1] =", self.points[1]
            print error

    def FindErrors(self):
        n = 10
        while True:
            self.Reset()
            self.GenerateData(n)

            try:
                self.Run()
            except Exception as error:
                print "!!! ERROR !!!"
                print "self.points[0] =", self.points[0]
                print "self.points[1] =", self.points[1]
                print error
                exit(1)

    def ExperimentTime(self) :
        fout = open('res.txt', 'w')
        nPoints = [100 * x for x in xrange(1, 31)]
        counter = 0
        numIter = 1
        for n in nPoints :
            while True :
                counter = counter + 1
                print "iter", counter
                self.Reset()
                self.GenerateData(n)

                try :
                    fictiveTime = self.Run()
                    break
                except Exception as error:
                    print error
                    pass
            fout.write(str(len(self.points[2])) + ' ' + str(self.scipyTime / numIter) + ' ' + str(fictiveTime / numIter) + '\n')
