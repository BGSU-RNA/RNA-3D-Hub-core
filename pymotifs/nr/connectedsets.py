# connectedsets.py takes as input a dictionary of sets.  The keys of the dictionary and the elements of the lists are the "vertices" and the presence of a vertex in a set means that the dictionary key and the vertex are linked
# connectedsets.py finds the connected sets
# it returns a dictionary indexed by the first structure found in each set.  There is nothing special about the keys of this dictionary.
# Note that connections[i] does not need to contain i; that is assumed
# Note that connections[i] may contain j without connections[j] containing i; the program adds i to connections[j]; it assumes symmetry.

def findconnectedsets(connections):

  considered = {}
  for i in connections.keys():
    considered[i] = False
    for j in connections[i]:
        considered[j] = False

  for i in connections.keys():
    connections[i] = connections[i] | set([i])    # make sure all connections are reflexive; i is connected to i
    for j in connections[i]:
        if j in connections:
            connections[j] = connections[j] | set([i])  # make sure all connections are entered in reversed order too; symmetrize
        else:
            connections[j] = set([i])

  linked = {}

  for i in connections.keys():              # loop through all keys for connections
    if not considered[i]:                   # if this key was not already encountered,
        linked[i] = connections[i]            # start a new list of its links
        considered[i] = True                  # note that it has been encountered, so it is not processed again
        newconnections = True                 # set a flag for exploring what i is connected to
        while newconnections:                 # keep doing this as long as a new connection was found
            newconnections = False              # no new connections found yet
            for j in linked[i]:                 # look through the things connected to i
                if not considered[j]:             # that have not already been considered
                    linked[i] = linked[i] | connections[j] # add these to the connections
                    considered[j] = True            # don't consider these connections again
                    newconnections = True           # a new connection was found
        print linked[i]

    return linked

if __name__ == "__main__":
    connections = {}
    connections['A'] = ['B','C']
    connections['B'] = ['D']
    connections['C'] = ['E']
    connections['E'] = ['A']
    connections['F'] = ['B']
    connections['zA'] = ['zB','zC']
    connections['zB'] = ['zD']
    connections['zC'] = ['zE']
    connections['zD'] = ['zF']
    connections['zE'] = ['zA']
    connections['zF'] = ['zB']
    for key in connections.keys():
        connections[key] = set(connections[key])


    findconnectedsets(connections)
