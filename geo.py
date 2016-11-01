import math

EPS = 1e-4
PI = math.pi
MINDIST = 2

def Sign(x):
    return int(math.copysign(1, x))

def GetLength(p1, p2):
    return ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2) ** (0.5)

def ExteriorProd(p1, p2, p3):
    v1 = [p1[0] - p2[0], p1[1] - p2[1]]
    v2 = [p3[0] - p2[0], p3[1] - p2[1]]

    res = v1[0] * v2[1] - v2[0] * v1[1]
    return res

def CountAngle(p1, p2, p3):
    l1 = GetLength(p1, p2)
    l2 = GetLength(p2, p3)
    
    if l1 * l2 < EPS:
        return 0

    v1 = [p2[0] - p1[0], p2[1] - p1[1]]
    v2 = [p2[0] - p3[0], p2[1] - p3[1]]

    sin_arg = (v1[0] * v2[1] - v2[0] * v1[1]) / (l1 * l2)
    if abs(abs(sin_arg) - 1) < EPS:
        x = Sign(sin_arg)
        aSin = math.asin(x)
    else:
        aSin = math.asin(sin_arg)

    cos_arg = (v1[0] * v2[0] + v1[1] * v2[1]) / (l1 * l2)
    if abs(abs(cos_arg) - 1) < EPS:
        x = Sign(cos_arg)
        aCos = math.acos(x)
    else :
        aCos = math.acos(cos_arg)
    
    if sin_arg >= 0 and cos_arg >= 0: # first quarter
        return (aSin + aCos) / 2.0
    elif sin_arg >= 0 and cos_arg < 0: # second quarter
        return ((PI - aSin) + aCos) / 2.0
    elif sin_arg < 0 and cos_arg < 0: # third quarter
        return ((PI - aSin) + (2 * PI - aCos)) / 2.0
    else: # fourth quarter
        return 2 * PI + (aSin - aCos) / 2.0

def LineParameters(p1, p2):
    a = p1[1] - p2[1]
    b = p2[0] - p1[0]
    c = p1[0] * p2[1] - p2[0] * p1[1]
    return [a, b, c]

def PerpendicularParameters(p, line):
    [a, b, c] = line
    pa = b
    pb = -a
    pc = -b * p[0] + a * p[1]
    return [pa, pb, pc]

def IntersectLines(line1, line2):
    [a1, b1, c1] = line1
    [a2, b2, c2] = line2
    det = a1 * b2 - a2 * b1
    if abs(det) < EPS:
        return [1000, 1000]
    x0 = (c2 * b1 - c1 * b2) / det
    y0 = (c2 * a1 - c1 * a2) / (-det)
    return [x0, y0]

def RadiusByLineAndPoint(line, radius, center, openPoint, checkPoint):
    [a, b, c] = LineParameters(openPoint, checkPoint)

    midlePoint = [round((openPoint[0] + checkPoint[0]) / 2), round((openPoint[1] + checkPoint[1]) / 2)]
    [pa, pb, pc] = PerpendicularParameters(midlePoint, [a, b, c])

    pIntersect = IntersectLines(line, [pa, pb, pc])

    if GetLength(pIntersect, center) < radius :
        return GetLength(pIntersect, checkPoint)
    else :
        return -1