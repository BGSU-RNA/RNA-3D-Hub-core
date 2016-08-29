import os
import json
import collections
from copy import deepcopy


def merge(a, b):
    """Recursively merges dict's. not just simple a['key'] = b['key'], if
    both a and b have a key who's value is a dict then merge is
    called on both values and the result stored in the returned
    dictionary.

    Taken from:
    https://www.xormedia.com/recursively-merge-dictionaries-in-python/

    :a: A dictionary.
    :b: Another dictionary.
    """

    if not isinstance(b, dict):
        return b

    result = deepcopy(a)
    for k, v in b.iteritems():
        if k in result and isinstance(result[k], dict):
            result[k] = merge(result[k], v)
        else:
            new_v = deepcopy(v)
            if isinstance(new_v, basestring):
                new_v = str(new_v)
            result[k] = new_v

    return result


def defaults():
    """A function to create a dictionary with the default values for the
    pipeline.
    """

    here = os.path.dirname(__file__)
    base = os.path.abspath(os.path.join(here, '..'))

    return {
        "locations": {
            "base": base,
            "cache": os.path.join(base, "cache"),
            "loops_mat_files": os.path.join(base, "MotifAtlas",
                                            "PrecomputedData"),
            "loops_search_dir": os.path.join(base, "MotifAtlas", "aAa"),
            "log_dir": os.path.join(base, "MotifAtlas", "logs"),
            "releases_dir": os.path.join(base, "MotifAtlas", "Releases"),
            "fr3d_root": os.path.join(base, "FR3D"),
            "2ds_destination":  os.path.join(base, "MotifAtlas", "2d"),
        },
        'recaculate': collections.defaultdict(lambda: False)
    }


def load(filename):
    config = collections.defaultdict(dict)
    with open(filename, 'rb') as raw:
        config.update(merge(defaults(), json.load(raw)))
        return config
