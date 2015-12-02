def total_ordering(cls):
    """Class decorator that fills in missing ordering methods. For details on
    what this does see:
    https://docs.python.org/2.7/library/functools.html#functools.total_ordering.
    This is backported to work with Python 2.6. Classes using this must
    implement one of __lt__, __le__, __gt__, or __ge__. They should also
    implement a __eq__ and __ne__ method for completeness.

    Code from: http://code.activestate.com/recipes/576685/

    This version is taken from: https://github.com/kvesteri/total-ordering
    """

    convert = {
        '__lt__': [
            (
                '__gt__',
                lambda self, other: not (self < other or self == other)
            ),
            (
                '__le__',
                lambda self, other: self < other or self == other
            ),
            (
                '__ge__',
                lambda self, other: not self < other
            )],
        '__le__': [
            (
                '__ge__',
                lambda self, other: not self <= other or self == other
            ),
            (
                '__lt__',
                lambda self, other: self <= other and not self == other
            ),
            (
                '__gt__',
                lambda self, other: not self <= other
            )],
        '__gt__': [
            (
                '__lt__',
                lambda self, other: not (self > other or self == other)
            ),
            (
                '__ge__',
                lambda self, other: self > other or self == other
            ),
            (
                '__le__',
                lambda self, other: not self > other
            )],
        '__ge__': [
            (
                '__le__',
                lambda self, other: (not self >= other) or self == other
            ),
            (
                '__gt__',
                lambda self, other: self >= other and not self == other
            ),
            (
                '__lt__',
                lambda self, other: not self >= other
            )]
    }

    roots = set(dir(cls)) & set(convert)
    if not roots:
        raise ValueError('Must define at least one ordering operation')

    # prefer __lt__ to __le__ to __gt__ to __ge__
    root = max(roots)
    for opname, opfunc in convert[root]:
        if opname not in roots:
            opfunc.__name__ = opname
            opfunc.__doc__ = getattr(int, opname).__doc__
            setattr(cls, opname, opfunc)
    return cls
