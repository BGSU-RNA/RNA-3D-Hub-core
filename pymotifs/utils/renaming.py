"""This module contains some utitlities useful when renaming and altering data.
It provides a relatively straightforward way to rename and alter data from one
form to another. The data structures are a bit basic right now, but could be
improved easily enough.
"""

import collections as coll


class Wrapper(coll.namedtuple('Wrapper', ['initial', 'final', 'fn'])):
    def __call__(self, *args, **kwargs):
        return self.fn(*args, **kwargs)


def maybe_str(value, strip=False, also_none=set()):
    val = value
    if val:
        if strip and isinstance(val, str):
            val = value.strip()
        if val and val.lower() not in also_none:
            return str(val)
    return None


def none_or(func, strip=True, also_none=set()):
    def fn(value, strip=strip, also_none=also_none):
        val = maybe_str(value, strip=strip, also_none=also_none)
        if not val:
            return None
        return func(val)
    return fn

maybe_float = none_or(float)
maybe_int = none_or(int)


def rename(name, func, **kwargs):
    def fn(data):
        return func(data.get(name, None), **kwargs)
    return Wrapper(initial=name, final=None, fn=fn)


def transform(name, func, **kwargs):
    return rename(name, func, **kwargs)._replace(final=name)


def with_dashes(name, func, **kwargs):
    pretty = name.replace('-', '_').lower()
    return rename(name, func, **kwargs)._replace(final=pretty)


class Renamer(object):
    def __init__(self, *args, **kwargs):
        self.renamers = []
        self.renamers.extend(args)

        for key, value in kwargs.items():
            wrapper = Wrapper(initial=value.initial, final=key, fn=value.fn)
            self.renamers.append(wrapper)

    def __call__(self, data, skip_missing=False):
        result = {}
        for wrapper in self.renamers:
            if wrapper.initial not in data and skip_missing:
                continue
            result[wrapper.final] = wrapper.fn(data)
        return result
