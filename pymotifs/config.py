"""A module for loading configuration data. This has a function to load and
apply the default values for pipeline configuration. It does not check if the
configuration is valid however.
"""

import os
import json
import collections
from copy import deepcopy
from six import iteritems, string_types


def merge(a, b):
    """Recursively merges dict's. not just simple a['key'] = b['key'], if
    both a and b have a key who's value is a dict then merge is
    called on both values and the result stored in the returned
    dictionary. This will not modify either input.

    Taken from:
    https://www.xormedia.com/recursively-merge-dictionaries-in-python/

    Parameters
    ----------
    a : dict
        A dictionary.
    b : dict
        Another dictionary.

    Returns
    -------
    merged : dict
        The merge dictonary.
    """

    if not isinstance(b, dict):
        return b

    result = deepcopy(a)
    for k, v in iteritems(b):
        if k in result and isinstance(result[k], dict):
            result[k] = merge(result[k], v)
        else:
            new_v = deepcopy(v)
            if isinstance(new_v, string_types):
                new_v = str(new_v)
            result[k] = new_v

    return result


def defaults():
    """A function to create a dictionary with the default values for the
    pipeline. Notably, this will compute the paths to most places as well as
    setting recaculate for all stages to False.

    Returns
    -------
    defaults : dict
        A dictonary of default values.
    """

    here = os.path.dirname(__file__)
    base = os.path.abspath(os.path.join(here, '..'))

    return {
        'db': {
            'pool_size': 40,
            'max_overflow': 10,
        },
        "locations": {
            "base": base,
            "cache": os.path.join(base, "cache"),
            "loops_mat_files": os.path.join(base, "MotifAtlas",
                                            "PrecomputedData"),
            "loops_search_dir": os.path.join(base, "MotifAtlas", "aAa"),
            "log_dir": os.path.join(base, "MotifAtlas", "logs"),
            "releases_dir": os.path.join(base, "MotifAtlas", "Releases"),
            "fr3d_root": os.path.join(base, "FR3D"),      # where .cif files are saved
            'quality_reports': os.path.join(base, "MotifAtlas", 'quality',
                                            'validation-reports'),
        },
        'recaculate': collections.defaultdict(lambda: False)
    }


def load(filename):
    """Load the configuration from the given file. This will load the JSON file
    at the specified path as well as merge in the default values.

    Parameters
    ----------
    filename : str
        The filename to load configuration from.

    Returns
    -------
    configuration : collections.defaultdict
        The configuration.
    """

    config = collections.defaultdict(dict)
    with open(filename, 'rb') as raw:
        config.update(merge(defaults(), json.load(raw)))
        return config
