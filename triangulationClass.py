from Tkinter import *
from math import *
from random import *
from time import *
from scipy.spatial import Delaunay

class triangulationClass() :
    root = Tk()
    
    readMode = True
    points = []
    faces = []
    edges = []
    openEdges = []
            
    triangleDone = False
    experimentMode = False

    def __init__(self, _width = 500, _height = 300) :
        self.width = _width
        self.height = _height
        self.root.title("Construction of triangulation")
        self.root.minsize(width = _width, height = _height)
        self.root.maxsize(width = _width, height = _height)

        self.canvas = Canvas(self.root, width = _width, height = _height, bg = "white")       
        self.canvas.pack()

        self.root.bind('<Button-1>', self.clickMouse)

        self.b1 = Button(self.root, bg = "white", fg = "blue", text = "O(n^4)", command = self.makeTriangle_n4)
        self.b1.place(x = 50, y = _height - 50)

        self.b2 = Button(self.root, bg = "white", fg = "blue", text = "O(n^3)", command = self.makeTriangle_n3)
        self.b2.place(x = 150, y = _height - 50)

        self.b3 = Button(self.root, bg = "white", fg = "blue", text = "O(n^2)", command = self.makeTriangle_n2)
        self.b3.place(x = 250, y = _height - 50)

        self.b4 = Button(self.root, bg = "white", fg = "blue", text = "O(nlogn)", command = self.makeTriangle_nlogn)
        self.b4.place(x = 350, y = _height - 50)

        self.b5 = Button(self.root, bg = "white", fg = "blue", text = "Experiment", command = self.experiment)
        self.b5.place(x = 450, y = _height - 50)

    def clickMouse(self, event) :
        if self.readMode == True :
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
        circle = self.canvas.create_oval(x - rad, y - rad, x + rad, y + rad, fill = "black")

    def drawTriangle(self) :
        if self.experimentMode == True :
            return
        
        for i in range(len(self.faces)) :
            temp = [self.points[self.faces[i][0]], self.points[self.faces[i][1]], self.points[self.faces[i][2]], self.points[self.faces[i][0]]]
            for j in range(3) :
                line = self.canvas.create_line(temp[j][0], temp[j][1], temp[j + 1][0], temp[j + 1][1])
    
    def checkFace(self, x) :
        x = sorted(x)
        return x not in self.faces
        
    def addFace(self, x) :
        x = sorted(x)
        if x not in self.faces :
            self.faces.append([x[0], x[1], x[2]])

    def addEdge(self, x) :
        temp = [x[1], x[0]]
        if temp not in self.edges :
            self.edges.append(x)
            self.openEdges.append(True)
        else :
            i = self.edges.index(temp)
            self.openEdges[i] = False
    
    def makeTriangle_n4(self) :
        if self.triangleDone == False :
            self.triangleDone = True    
            if self.experimentMode == False :
                self.points.pop(len(self.points) - 1)
            self.readMode = False

            numPoints = len(self.points)

            for i in range(0, numPoints - 2) :
                for j in range(i + 1, numPoints - 1) :
                    for k in range(j + 1, numPoints) :
                        res = self.checkDeloneConditionCircle(i, j, k)
                        if res == True :
                            self.addFace([i, j, k])

            self.drawTriangle()
    
    def checkDeloneConditionCircle(self, n1, n2, n3) :
        res = self.centerCircle(self.points[n1], self.points[n2], self.points[n3])
        if res == False :
            return False
    
        center = res[0]
        rad = res[1]
    
        flag = True
        for i in range(0, len(self.points)) :
            if (i != n1 and i != n2 and i != n3) :
                if (center[0] - self.points[i][0]) ** 2 + (center[1] - self.points[i][1]) ** 2 < rad ** 2 :
                    #return False
                    flag = False
        return flag
    
    def centerCircle(self, p1, p2, p3) :
        a = p2[0] - p1[0]
        b = p2[1] - p1[1]
        c = p3[0] - p1[0]
        d = p3[1] - p1[1]

        e = a * (p1[0] + p2[0]) + b * (p1[1] + p2[1])
        f = c * (p1[0] + p3[0]) + d * (p1[1] + p3[1])
        g = 2 * (a * (p3[1] - p2[1]) - b * (p3[0] - p2[0]))

        if g == 0 :
            return False

        centerX = (d * e - b * f) / g
        centerY = (a * f - c * e) / g
        center = [centerX, centerY]
        radius = ((center[0] - p1[0]) ** 2 + (center[1] - p1[1]) ** 2) ** 0.5

        return [center, radius]

    def makeTriangle_n3(self) :
        if self.triangleDone == False :
            self.triangleDone = True    
            if self.experimentMode == False :
                self.points.pop(len(self.points) - 1)
            self.readMode = False

            numPoints = len(self.points)

            for i in range(0, numPoints - 1) :
                for j in range(i + 1, numPoints) :
                    self.checkDeloneConditionLine(i, j)

            self.drawTriangle()
    
    def checkDeloneConditionLine(self, n1, n2) :
        [a, b, c] = self.lineParameters(self.points[n1], self.points[n2])
        angles = [False, False] # 0 - sign < 0, 1 - sign > 0
        iAngles = [-1, -1]
        for i in range(0, len(self.points)) :
            if (i != n1 and i != n2) :
                sign = a * self.points[i][0] + b * self.points[i][1] + c
                temp = self.countAngle(self.points[n1], self.points[i], self.points[n2])

                if sign == 0 :
                    continue
                elif sign > 0 :
                    sign = 1
                else :
                    sign = 0
                if angles[sign] == False or angles[sign] < temp : 
                    angles[sign] = temp
                    iAngles[sign] = i
        
        if angles[0] == False and angles[1] == False :
            return
        elif angles[0] == False and angles[1] != False :
            self.addFace([n1, n2, iAngles[1]])
        elif angles[0] != False and angles[1] == False :
            self.addFace([n1, n2, iAngles[0]])
        elif angles[0] + angles[1] < acos(-1.0) :
            self.addFace([n1, n2, iAngles[0]])
            self.addFace([n1, n2, iAngles[1]])
        
    def lineParameters(self, p1, p2) :
        a = p1[1] - p2[1]
        b = p2[0] - p1[0]
        c = p1[0] * p2[1] - p2[0] * p1[1]
        return [a, b, c]

    def countAngle(self, p1, p2, p3) :
        v1 = [p2[0] - p1[0], p2[1] - p1[1]]
        v2 = [p2[0] - p3[0], p2[1] - p3[1]]

        res = acos((v1[0] * v2[0] + v1[1] * v2[1]) / ((v1[0] ** 2 + v1[1] ** 2) * (v2[0] ** 2 + v2[1] ** 2)) ** 0.5)
        return res

    def makeTriangle_n2(self) :
        if self.triangleDone == False :
            self.triangleDone = True    
            if self.experimentMode == False :
                self.points.pop(len(self.points) - 1)
            self.readMode = False

            numPoints = len(self.points)
            
            numBottom = 0
            for i in range(1, numPoints) :
                if self.points[i][1] > self.points[numBottom][1] :
                    numBottom = i
                elif self.points[i][1] == self.points[numBottom][1] and self.points[i][0] > self.points[numBottom][0] :
                    numBottom = i
            
            tempPoint = [self.points[numBottom][0] + 5, self.points[numBottom][1]]
            angle = False
            iAngle = -1                    
            for i in range(0, numPoints) :
                if i != numBottom :
                    temp = self.countAngle(self.points[i], self.points[numBottom], tempPoint)
                    if angle == False or angle > temp:
                        angle = temp
                        iAngle = i
            
            self.addEdge([numBottom, iAngle])
            
            i = -1
            while(True) :
                i += 1
                if i == len(self.edges) :
                    break
                if self.openEdges[i] == True :
                    self.openEdges[i] = False
                    n1 = self.edges[i][0]
                    n2 = self.edges[i][1]
                    
                    angle = False
                    iAngle = -1                    
                    for j in range(0, len(self.points)) :
                        if j == n1 or j == n2 :
                            continue

                        v1 = [self.points[n2][0] - self.points[n1][0], self.points[n2][1] - self.points[n1][1]]
                        v2 = [self.points[j][0] - self.points[n1][0], self.points[j][1] - self.points[n1][1]]
                        
                        sign = v1[0] * v2[1] - v1[1] * v2[0]
                        if self.checkFace([n1, n2, j]) and sign < 0 :
                            temp = self.countAngle(self.points[n1], self.points[j], self.points[n2])
                            if angle == False or angle < temp :
                                angle = temp
                                iAngle = j

                    if angle :
                        self.addFace([n1, n2, iAngle])
                        self.addEdge([n1, iAngle])
                        self.addEdge([iAngle, n2])
                      
            self.drawTriangle()
    
    def makeTriangle_nlogn(self) :
        if self.triangleDone == False :
            self.triangleDone = True    
            if self.experimentMode == False :
                self.points.pop(len(self.points) - 1)
            self.readMode = False

            tri = Delaunay(self.points)
            
            for f in tri.vertices :
                self.faces.append([f[0], f[1], f[2]])
        
            self.drawTriangle()

    """def makeTriangle_nlogn(self) :
        if self.triangleDone == False :
            self.triangleDone = True    
            if self.experimentMode == False :
                self.points.pop(len(self.points) - 1)
            self.readMode = False

            
            self.drawTriangle()"""
                
    def experiment(self) :
        if self.experimentMode == False :
            self.experimentMode = True
            self.readMode = False

            fout = open('res.txt', 'w')

            nPoints = [20, 30, 50, 70, 90, 110]
            modes = ["O(nlogn)", "O(n^2)", "O(n^3)", "O(n^4)"]
            result = []

            for n in nPoints :
                fout.write(str(n) + ' ')
                self.points = []
                for i in range(n) :
                    x = randrange(self.width - 10) + 5
                    y = randrange(self.height - 60) + 5
 
                    self.points.append([x, y])
                
                temp = []
                for i in modes :
                    res = self.make_experiment(i)
                    fout.write(str(res) + ' ')
                    temp.append(res)

                result.append(temp)
                fout.write('\n')

            self.experimentMode = False
            self.drawTriangle()
            self.experimentMode = True
            print("Done")

    def make_experiment(self, mode) :
        self.faces = []
        self.edges = []
        self.openEdges = []
        self.triangleDone = False
        
        start = time()
        
        if mode == "O(n^4)" :
            self.makeTriangle_n4()
        elif mode == "O(n^3)" :
            self.makeTriangle_n3()
        elif mode == "O(n^2)" :
            self.makeTriangle_n2()
        elif mode == "O(nlogn)" :
            self.makeTriangle_nlogn()

        finish = time()

        print(mode, 'points', len(self.points), 'time', finish - start)
        return finish - start