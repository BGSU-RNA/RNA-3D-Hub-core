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


class BadlyFormattedRelease(Exception):
    """This is raised if we are given a release id that is not formatted
    correctly. A correctly formatted id must have two parts, major and minor
    where each is a integer.
    """
    pass


class ImpossibleRelease(Exception):
    """This is raised if we are asking for a preivous release which is
    impossible. For example release -1.0 would be an example which can happen
    when getting the previous release of 0.0.
    """
    pass


class Release(object):

    names = ['major', 'minor']

    def __init__(self, config, session):
        self.session = core.Session(session)
        self.config = config

    def previous(self, release, mode='minor', bad_id='wrap'):
        if mode not in self.names:
            raise UnknownReleaseMode(mode)

        if bad_id not in ['wrap', 'none', 'raise']:
            raise ValueError('Invalid bad_id action %s' % bad_id)

        current = self.__parse__(release)
        if mode == 'minor':
            current['minor'] -= 1

        if mode == 'major':
            current['minor'] = 0
            current['major'] -= 1

        if current['major'] < 0:
            if bad_id == 'none':
                return None

            if bad_id == 'wrap':
                current['major'] = 0
                current['minor'] = 1

            if bad_id == 'raise':
                raise ImpossibleRelease(release)

        if current['minor'] < 0:
            if bad_id == 'wrap':
                current['minor'] = 0

            if bad_id == 'none':
                return None

            if bad_id == 'raise':
                raise ImpossibleRelease(release)

        if current['major'] == 0 and current['minor'] == 0:
            current['minor'] = 1

        return self.__format__(current)

    def current(self, name):
        if name not in ['nr', 'loop', 'motif']:
            raise UnknownReleaseType(name)

        tablename = 'current_%s_release' % name
        with self.session() as session:
            result = session.execute(CURRENT_QUERY % tablename)
            result = result.fetchall()
            if len(result) == 0:
                return '0.0'
            return result[0].id

    def next(self, current, mode='minor'):
        if mode == 'none' or mode is None:
            return str(current)

        if mode not in self.names:
            raise UnknownReleaseMode(mode)

        current = self.__parse__(current)
        if mode == 'minor':
            current['minor'] += 1

        if mode == 'major':
            current['minor'] = 0
            current['major'] += 1

        return self.__format__(current)

    def __parse__(self, release_id):
        try:
            parts = map(int, release_id.split('.'))
        except:
            raise BadlyFormattedRelease("Can't process release id: %s" %
                                        release_id)

        if len(parts) != 2:
            raise BadlyFormattedRelease("Release id must have major an minor")

        return dict(zip(self.names, parts))

    def __format__(self, current):
        return '.'.join((str(current['major']), str(current['minor'])))
