#from scipy.spatial import Delaunay
#points = [[171, 56], [295, 57], [223, 94], [95, 89], [179, 110], [349, 116], [293, 132], [60, 203], [124, 212], [165, 167], [247, 172], [216, 231], [332, 201]]

#tri = Delaunay(points)

#print(tri.simplices)

#for i in tri.simplices :
#    print(i)

#self.points = [[171, 56], [295, 57], [223, 94], [95, 89], [179, 110], [349, 116], [293, 132], [60, 203], [124, 212], [165, 167], [247, 172], [216, 231], [332, 201], [1, 1]]
 
 # coding: utf-8
from Tkinter import *

import time
def deltext():
    canvas.delete(text)
def delsq():
    canvas.delete(square)
    
root = Tk()
btn1 = Button(root, text="Delete text", command=deltext)
btn1.pack()
btn2 = Button(root, text="Delete square", command=delsq)
btn2.pack()
canvas = Canvas(root, width=300, height=300)
canvas.pack(fill=BOTH)
circle = canvas.create_oval(0,0,300,300, fill="blue")
diamond = canvas.create_polygon(150,0,0,150,150,300,300,150, fill="red")
square = canvas.create_rectangle(75,75,225,225, fill="green")

text = canvas.create_text(150,150, text="Tkinter canvas", fill="purple", 
        font=("Helvectica", "16"))
root.mainloop()