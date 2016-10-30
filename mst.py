from geo import getLength

def ConstructMST(points, faces):
    edgesMST = [[], []]
    tree = [[], []]
        
    for i in xrange(2):
        for f in faces[i] :
            p = [f[0], f[1], f[2], f[0]]
            
            for j in xrange(3):
                lenEdge = getLength(points[i][p[j]], points[i][p[j + 1]])

                if edgesMST[i].count((lenEdge, p[j], p[j + 1])) == 0 :
                    edgesMST[i].append((lenEdge, p[j], p[j + 1]))

        edgesMST[i].sort()
            
        treeId = []
        for j in xrange(len(points[i])):
            treeId.append(j)

        for e in edgesMST[i]:
            a = e[1]
            b = e[2]

            if treeId[a] != treeId[b]:
                tree[i].append((a, b))
                oldId = treeId[b]
                newId = treeId[a]
                for j in xrange(len(treeId)):
                    if treeId[j] == oldId:
                        treeId[j] = newId

            
    mst = [[], [], []]
    n1 = len(points[0])
        
    for i in xrange(2):
        for j in xrange(len(tree[i])):
            mst[i].append([tree[i][j][0] + i * n1, tree[i][j][1] + i * n1])

    mst[2].extend(mst[0])
    mst[2].extend(mst[1])
    return mst