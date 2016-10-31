from Tkinter import Tk, Canvas, Button

import draw

class BaseForm(object):
    # input first set of points
    firstPointsSet = True

    preprocessingDone = False

    # in experement mode do not remove last point (click on button)
    experimentMode = False

    points = [[], [], []]

    # data for creation structure
    faces = [[], [], []]
    neighFaces = [[], [], []]
    
    # data for detecting next starter
    fictiveEdges = []
    
    scipyTime = -1

    def __init__(self, _width = 600, _height = 400):
        self.root = Tk()
        self.width = _width
        self.height = _height
        self.root.title("Construction of triangulation")
        self.root.minsize(width = _width, height = _height)
        self.root.maxsize(width = _width, height = _height)

        self.canvas = Canvas(self.root, width = _width, height = _height, bg = "white")       
        self.canvas.pack()

        self.root.bind('<Button-1>', self.ClickMouse)

        self.b1 = Button(self.root, bg = "white", fg = "black", text = "Second points set", command = self.SecondPointsSet)
        self.b1.place(x = 50, y = _height - 50)        
        
        self.b2 = Button(self.root, bg = "white", fg = "black", text = "Triangulation", command = self.Run)
        self.b2.place(x = 200, y = _height - 50)
        
        self.b3 = Button(self.root, bg = "white", fg = "black", text = "Experiment", command = self.Experiment)
        self.b3.place(x = 325, y = _height - 50)
        
        self.b4 = Button(self.root, bg = "white", fg = "black", text = "Errors", command = self.FindErrors)
        self.b4.place(x = 425, y = _height - 50)
    
        self.b5 = Button(self.root, bg = "white", fg = "black", text = "Time", command = self.ExperimentTime)
        self.b5.place(x = 500, y = _height - 50)

    def ClickMouse(self, event):
        if not self.preprocessingDone:
            getX = event.x_root
            getY = event.y_root

            posRootX = self.root.winfo_rootx()
            posRootY = self.root.winfo_rooty()
    
            x = getX - posRootX
            y = getY - posRootY

            p = [x, y]
            self.points[not self.firstPointsSet].append(p)
            draw.Point(self.canvas, p, "black", 3, not self.firstPointsSet)
    
    def Reset(self):
        self.canvas.delete("all")

        self.firstPointsSet = True
        self.preprocessingDone = False
        self.experimentMode = True
        
        self.points = [[], [], []]

        self.faces = [[], [], []]
        self.neighFaces = [[], [], []]

        self.fictiveEdges = []
        
        self.scipyTime = -1

    def SecondPointsSet(self):
        self.firstPointsSet = False

    def Run(self):
        raise NotImplementedError('This function must be implemented within child class!')

    def Experiment(self):
        raise NotImplementedError('This function must be implemented within child class!')

    def FindErrors(self):
        raise NotImplementedError('This function must be implemented within child class!')

    def ExperimentTime(self):
        raise NotImplementedError('This function must be implemented within child class!')