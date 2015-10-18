from pymotifs import core
from pymotifs.models import MlReleases
from pymotifs.models import LoopReleases
from pymotifs.models import NrReleases

from sqlalchemy import desc


CURRENT_QUERY = """
SELECT {column}
FROM {table}
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


class Release(core.Base):

    names = ['major', 'minor', 'lookup']

    column_mapping = {
        'motif': 'ml_releases_id',
        'loop': 'loop_releases_id',
        'nr': 'nr_release_id'
    }

    table_mapping = {
        'motif': MlReleases,
        'loop': LoopReleases,
        'nr': NrReleases
    }

    def previous(self, release, mode='minor', bad_id='wrap', name=None):
        if mode not in self.names:
            raise UnknownReleaseMode(mode)

        if bad_id not in ['wrap', 'none', 'raise']:
            raise ValueError('Invalid bad_id action %s' % bad_id)

        current = self.__parse__(release)
        if mode == 'minor':
            current['minor'] -= 1

        if mode == 'major':
            current['major'] -= 1
            current['minor'] = 0

        if mode == 'lookup':
            column_name = self.column_mapping[name]
            table = self.table_mapping[name]
            column = getattr(table, column_name)
            with self.session() as session:
                query = session.query(column).\
                    order_by(desc(table.date)).\
                    limit(2)

                results = [getattr(result, column_name) for result in query]
                try:
                    index = results.index(release)
                    return results[index + 1]
                except:
                    return None

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

        column_name = self.column_mapping[name]
        table = self.table_mapping[name]
        column = getattr(table, column_name)
        with self.session() as session:
            query = session.query(column).\
                order_by(desc(table.date)).\
                limit(1)

            if query.count() == 0:
                return '0.0'

            return getattr(query.one(), column_name)

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
