from math import *

EPS = 1e-4
PI = pi

def getLength(p1, p2) :
    return ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2) ** (0.5)

def countSignAngle(p1, p2, p3) :
    v1 = [p1[0] - p2[0], p1[1] - p2[1]]
    v2 = [p3[0] - p2[0], p3[1] - p2[1]]

    res = v1[0] * v2[1] - v2[0] * v1[1]
    return res

def countAngle(p1, p2, p3) :
    l1 = getLength(p1, p2)
    l2 = getLength(p2, p3)
        
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

def lineParameters(p1, p2) :
    a = p1[1] - p2[1]
    b = p2[0] - p1[0]
    c = p1[0] * p2[1] - p2[0] * p1[1]
    return [a, b, c]

def perpendicularParameters(p, line) :
    [a, b, c] = line
    pa = b
    pb = -a
    pc = -b * p[0] + a * p[1]
    return [pa, pb, pc]

def intersectLines(line1, line2) :
    [a1, b1, c1] = line1
    [a2, b2, c2] = line2
    det = a1 * b2 - a2 * b1
    if abs(det) < EPS :
        return [1000, 1000]
    x0 = (c2 * b1 - c1 * b2) / det
    y0 = (c2 * a1 - c1 * a2) / (-det)
    return [x0, y0]

def inSegment(p1, p2, checkPoint) :
    return min(p1[0], p2[0]) <= checkPoint[0] and checkPoint[0] <= max(p1[0], p2[0]) and \
        min(p1[1], p2[1]) <= checkPoint[1] and checkPoint[1] <= max(p1[1], p2[1])