#######################################################################
# Implements a topological sort algorithm.
#
# Copyright 2014 True Blade Systems, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Notes:
#  Modified from:
#     https://bitbucket.org/ericvsmith/toposort/
#  with some modifications for 2.6 comptability found in:
#      https://bitbucket.org/ericvsmith/toposort/pull-requests/1/fixes-for-python-26-compatibility
#  and some minor reformatting.
#  Based on http://code.activestate.com/recipes/578272-topological-sort
#   with these major changes:
#    Added unittests.
#    Deleted doctests (maybe not the best idea in the world, but it cleans
#     up the docstring).
#    Moved functools import to the top of the file.
#    Changed assert to a ValueError.
#    Changed iter[items|keys] to [items|keys], for python 3
#     compatibility. I don't think it matters for python 2 these are
#     now lists instead of iterables.
#    Copy the input so as to leave it unmodified.
#    Renamed function from toposort2 to toposort.
#    Handle empty input.
#    Switch tests to use set literals.
#
########################################################################

from functools import reduce as _reduce
from functools import partial

__all__ = ['levels', 'toposort']


def levels(data):
    """Dependencies are expressed as a dictionary whose keys are items
    and whose values are a set of dependent items. Output is a list of
    sets in topological order. The first set consists of items with no
    dependences, each subsequent set consists of items that depend upon
    items in the preceeding sets.
    """

    # Special case empty input.
    if len(data) == 0:
        return

    # Copy the input so as to leave it unmodified.
    data = data.copy()

    # Ignore self dependencies.
    for k, v in data.items():
        v.discard(k)

    # Find all items that don't depend on anything.
    extra_items_in_deps = _reduce(set.union, data.values()) - set(data.keys())

    # Add empty dependences where needed.
    for item in extra_items_in_deps:
        data[item] = set()

    while True:
        ordered = set(item for item, dep in data.items() if len(dep) == 0)
        if not ordered:
            break
        yield ordered
        data = dict(
            (item, (dep - ordered))
            for item, dep in data.items()
            if item not in ordered
        )

    if len(data) != 0:
        dat = ', '.join(repr(x) for x in data.items())
        msg = 'Cyclic dependencies exist among these items: {}'
        raise ValueError(msg.format(dat))


def toposort(data, by=None):
    """
    Returns a single list of dependencies. For any set returned by
    toposort(), those items are sorted and appended to the result (just to
    make the results deterministic).

    Not sure that this runs when you run the pipeline.
    """

    result = []
    for d in levels(data):
        fn = list
        if by:
            fn = partial(sorted, key=by)
        result.extend(fn(d))
    return result
