from Tkinter import *
from math import *
from random import *
from time import *
#from scipy.spatial import Delaunay

class triangulationClass() :
    root = Tk()
    
    readMode = True
    points = []
    
    def __init__(self, _width = 600, _height = 400) :
        self.width = _width
        self.height = _height
        self.root.title("Construction of triangulation")
        self.root.minsize(width = _width, height = _height)
        self.root.maxsize(width = _width, height = _height)

        self.canvas = Canvas(self.root, width = _width, height = _height, bg = "white")       
        self.canvas.pack()

        self.root.bind('<Button-1>', self.clickMouse)

        self.b1 = Button(self.root, bg = "white", fg = "blue", text = "Spanning tree", command = self.makeSpanningTree)
        self.b1.place(x = 50, y = _height - 50)

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

    def makeSpanningTree(setf) :
        pass