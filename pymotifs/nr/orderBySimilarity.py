import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import random
import math

def treePenalty(distance,link="average"):

    """
    print(distance.shape)

    for i in range(0,10):
        for j in range(0,10):
            print("%5.3f" % distance[i][j]),
        print
    """


    Z = linkage(squareform(distance),link)
#    print("regular",Z)

    penalty = np.zeros(distance.shape)

    group = []
    for i in range(0,distance.shape[0]):
        group.append([i])

    for merger in Z:
        a = int(merger[0])
        b = int(merger[1])
        group.append(group[a]+group[b])
        for i in group[a]:
            for j in group[b]:
                penalty[i][j] = merger[2]
                penalty[j][i] = merger[2]

    if 0 > 1:
        for i in range(0,len(distance)):
            for j in range(0,len(distance)):
                print("%6.3f" % penalty[i][j]),
            print("")


    # scale the penalty as appropriate.  Avoid dividing by 0.
    penalty = penalty * np.mean(distance) / max(0.00000001,np.mean(penalty))

    return penalty

def greedyInsertionPathLength(distance, order=[], verbose=False):
    # if no starting ordering
    if len(order) == 0:
        order = list(range(0,len(distance)))
        random.shuffle(order)          # random starting ordering
    path = order[:2]            # first two points of the current ordering
    score = distance[path[0], path[1]]

    for p in range(2, len(order)):
        # score inserting point order[p] at beginning of path
        bestScore = distance[order[p]][path[0]]
        bestPosition = 0

        # score inserting point order[p] at end of path
        currentScore = distance[path[-1]][order[p]]
        if currentScore < bestScore:
            bestScore = distance[path[-1]][order[p]]
            bestPosition = len(path)

        # score inserting point order[p] at various points within the path
        for position in range(1, len(path)):
            currentScore = distance[path[position-1]][order[p]] + distance[order[p]][path[position]] - distance[path[position-1]][path[position]]

            if currentScore < bestScore:
                bestScore = currentScore
                bestPosition = position

        path.insert(bestPosition, order[p])
        score += bestScore

    return path, score

def multipleGreedyInsertionPathLength(distance, repetitions=100):

    bestScore = float("inf")
    for rep in range(0,repetitions):
        path, score = greedyInsertionPathLength(distance)
        if score < bestScore:
            bestScore = score
            bestPath = path

    return bestPath

# orient the path to have lower discrepancies in upper left of matrix
def orientPath(distance,path):
    n = distance.shape[0]
    if n <= 2:
        return path
    else:
        if n == 3:
            m = 2
        elif n == 4:
            m = 3
        elif n == 5:
            m = 3
        elif n == 6:
            m = 4
        else:
            m = int(math.ceil(n/3))+1
        upperLeft = 0
        lowerRight = 0
        for i in range(0,m):
            for j in range(i+1,m):
                upperLeft  += distance[path[i]][path[j]]
                lowerRight += distance[path[n-i-1]][path[n-j-1]]

        if upperLeft > lowerRight:
            reversedPath = [path[n-i-1] for i in range(0,n)]
            return reversedPath
        else:
            return path

def treePenalizedPathLength(distance,repetitions=100,seed=None):
    if seed:
        random.seed(seed)

    n = distance.shape[0]
    if n > 2:
        penalizedMatrix = distance + treePenalty(distance)
        order = multipleGreedyInsertionPathLength(penalizedMatrix,repetitions)
        order = orientPath(distance,order)
    else:
        order = list(range(0,n))
    return order

def reorderSymmetricMatrix(distance, newOrder):

    newDistance = np.zeros(distance.shape)
    for i in range(0,newDistance.shape[0]):
        for j in range(i+1,newDistance.shape[1]):
           newDistance[i][j] = distance[newOrder[i]][newOrder[j]]
           newDistance[j][i] = newDistance[i][j]

    return newDistance

def reorderList(oldList,newOrder):
    newList = []
    for i in range(0,len(oldList)):
        newList.append(oldList[newOrder[i]])

    return newList

def imputeNANValues(distance):
    newDistance = np.zeros(distance.shape)

    maxVal = 0
    for i in range(0,newDistance.shape[0]):
        for j in range(i+1,newDistance.shape[1]):
            if not math.isnan(distance[i][j]):
                maxVal = max(maxVal,distance[i][j])

    for i in range(0,newDistance.shape[0]):
        for j in range(i+1,newDistance.shape[1]):
            if math.isnan(distance[i][j]) or distance[i][j] < 0:
                newDistance[i][j] = maxVal
                newDistance[j][i] = maxVal
            else:
                newDistance[i][j] = distance[i][j]
                newDistance[j][i] = distance[i][j]

        newDistance[i][i] = 0

    return newDistance

def optimalLeafOrder(distance):
    Z = linkage(squareform(distance), "average", optimal_ordering = True)
    dn = dendrogram(Z,no_plot=True)

    return dn['leaves']

def testPenaltyMatrix():
    data = [0,1,3,5,7,10]
    distance = np.zeros((len(data),len(data)))
    for i in range(0,len(data)):
        for j in range(0,len(data)):
            distance[i][j] = abs(data[i]-data[j])

    print("Distance matrix")
    print(distance)

    for linkage in ["single","complete","average","ward"]:
        print("Tree penalty for linkage %s" % linkage)
        tP = treePenalty(distance,linkage)
        print(tP)

def generateUniformDataset(n,d,seed=None):
    if seed:
        np.random.seed(seed)

    points = np.random.uniform(0,1,(n,d))

    distance = np.zeros((n,n))
    for i in range(0,n):
        for j in range(i+1,n):
            distance[i][j] = np.linalg.norm(points[i]-points[j])
            distance[j][i] = distance[i][j]

    return points, distance

def testOrdering():
    size = 150
    points, distance = generateUniformDataset(size,3,2276393)
    order = treePenalizedPathLength(distance,20,39873)
    print("Points:")
    print(points)
    print("Ordering:")
    print(order)

#testOrdering()
#testPenaltyMatrix()
