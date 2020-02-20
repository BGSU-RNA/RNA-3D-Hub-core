#!/usr/local/bin/python -O

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
from numpy import median
from numpy import isnan
from random import shuffle


def fixDistanceMatrix(distance,scanForNan=False):

  if scanForNan:
    values = []
    for i in range(0,len(distance)):
      for j in range(0,len(distance[0])):
        if not isnan(distance[i,j]) and distance[i,j] >= 0:
          values.append(distance[i,j])

    m = 0
    if len(values) > 0:
      m = median(values)
      for i in range(0,len(distance)):
        distance[i,i] = 0
        for j in range(0,len(distance[0])):
          distance[j,i] = distance[i,j]
          if isnan(distance[i,j]) or distance[i,j] == -1:
            distance[i,j] = m
            distance[j,i] = m
          if distance[i,j] < 0:
            distance[i,j] = 0
            distance[j,i] = 0

    return distance

def orderEquivalenceClassWithOLO(distance,scanForNan=False):

  distance = fixDistanceMatrix(distance,scanForNan)

#  Z = linkage(squareform(distance), "average", optimal_ordering = True) # not available on the server 8/8/2019
  Z = linkage(squareform(distance), "average")
  dn = dendrogram(Z,no_plot=True)

  bestOrder = dn['leaves']

  print("ordering",len(distance),bestOrder)

  return bestOrder

def greedyInsertionPathLength(distance, order=[], verbose=False):

    # if no starting ordering
    if len(order) == 0:
        order = range(0,len(distance))
        shuffle(order)          # random starting ordering
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

def orderEquivalenceClassWithPathLength(distance,scanForNan=False,repetitions=100):

  distance = fixDistanceMatrix(distance,scanForNan)

  bestOrder = multipleGreedyInsertionPathLength(distance, repetitions)

  print("ordering",len(distance),bestOrder)

  return bestOrder

