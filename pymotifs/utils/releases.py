from pymotifs import core


CURRENT_QUERY = """
SELECT id
FROM %s
LIMIT 1;
"""


class UnknownReleaseType(Exception):
    """This is raised when we are asked for the current release of an unknown
    type of release.
    """
    pass


class UnknownReleaseMode(Exception):
    """This is raised when the requested release mode is not known.
    """
    pass


class Release(object):

    names = ['major', 'minor']

    def __init__(self, config, session):
        self.session = core.Session(session)
        self.config = config

    def current(self, name):
        if name not in ['nr', 'loop', 'motif']:
            raise UnknownReleaseType(name)

        tablename = 'current_%s_release' % name
        with self.session() as session:
            result = session.execute(CURRENT_QUERY % tablename)
            result = result.fetchall()
            return result[0].id

    def next(self, current, mode='minor'):
        if mode == 'none' or mode is None:
            return str(current)

        if mode not in self.names:
            raise UnknownReleaseMode(mode)

        parts = map(int, current.split('.'))
        current = dict(zip(self.names, parts))
        func = getattr(self, '__%s__' % mode)

        return self.__format__(func(current))

    def __major__(self, current):
        return {
            'major': current['major'] + 1,
            'minor': 0
        }

    def __minor__(self, current):
        return {
            'major': current['major'],
            'minor': current['minor'] + 1
        }

    def __format__(self, current):
        return '.'.join((str(current['major']), str(current['minor'])))
