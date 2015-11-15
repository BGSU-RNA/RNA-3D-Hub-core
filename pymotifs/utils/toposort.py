import itertools as it


def toposort(graph, by=None):
    """Dependencies are expressed as a dictionary whose keys are items
    and whose values are a set of dependent items. Output is a list of
    sets in topological order. The first set consists of items with no
    dependences, each subsequent set consists of items that depend upon
    items in the preceeding sets.
    """

    levels_by_name = {}
    names_by_level = dict()

    def walk_depth_first(name):
        if name in levels_by_name:
            return levels_by_name[name]
        children = graph.get(name, None)
        level = 0
        if children:
            level = (1 + max(walk_depth_first(lname) for lname in children))
        levels_by_name[name] = level
        if level not in names_by_level:
            names_by_level[level] = set()
        names_by_level[level].add(name)
        return level

    for name in graph:
        walk_depth_first(name)

    if by:
        for level, stages in names_by_level.items():
            names_by_level[level] = sorted(stages, key=by)

    ordered_levels = (names_by_level.get(i, None) for i in it.count())
    ordered_sets = it.takewhile(lambda x: x is not None, ordered_levels)
    return it.chain.from_iterable(ordered_sets)
