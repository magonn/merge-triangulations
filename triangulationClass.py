from Tkinter import *
from math import *
from random import *
from time import *
from scipy.spatial import Delaunay

class triangulationClass() :
    root = Tk()
    
    trianglationDone = False
    points = []
    faces = []
    
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
                circle = self.canvas.create_oval(x - rad, y - rad, x + rad, y + rad, fill = "black", outline = "black"))
            
    def drawTriangle(self) :
        for i in range(len(self.faces)) :
            temp = [self.points[self.faces[i][0]], self.points[self.faces[i][1]], self.points[self.faces[i][2]], self.points[self.faces[i][0]]]
            for j in range(3) :
                self.line = self.canvas.create_line(temp[j][0], temp[j][1], temp[j + 1][0], temp[j + 1][1])

    def makeTriangle(self) :
        if self.trianglationDone == False :
            self.trianglationDone = True
            
            self.points.pop(len(self.points) - 1)
            tri = Delaunay(self.points)
            
            for f in tri.vertices :
                self.faces.append([f[0], f[1], f[2]])
        
            self.drawTriangle()

    def makeMST(self) :
        self.makeTriangle()