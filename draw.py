def drawPoint(canvas, p, col, wid = 2, type = 0) :
    if type == 0 :
        circle = canvas.create_oval(p[0] - wid, p[1] - wid, p[0] + wid, p[1] + wid, fill = col, outline = col)
    else :
        circle = canvas.create_oval(p[0] - wid, p[1] - wid, p[0] + wid, p[1] + wid, fill = "white", outline = col)

def drawCircle(canvas, p, col, wid = 0, type = 0) :
    if type == 0 :
        circle = canvas.create_oval(p[0] - wid, p[1] - wid, p[0] + wid, p[1] + wid, outline = col)
    else :
        circle = canvas.create_oval(p[0] - wid, p[1] - wid, p[0] + wid, p[1] + wid, outline = col, dash = '- -')

def drawEdge(canvas, p1, p2, col, wid = 1, type = 0) :
    if type == 0 or type == 2 :
        line = canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill = col, width = wid)
    else :
        line = canvas.create_line(p1[0], p1[1], p2[0], p2[1], fill = col, width = wid, dash = '- -')
            
def drawTriangle(canvas, points, faces, col, wid = 1, type = 0) :
    for f in faces :
        temp = [points[f[0]], points[f[1]], points[f[2]], points[f[0]]]
            
        for j in range(3) :
            if type == 0 or type == 2:
                line = canvas.create_line(temp[j][0], temp[j][1], temp[j + 1][0], temp[j + 1][1], fill = col, width = wid)
            else :
                line = canvas.create_line(temp[j][0], temp[j][1], temp[j + 1][0], temp[j + 1][1], fill = col, width = wid, dash = '- -')
