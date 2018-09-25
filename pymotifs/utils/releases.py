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


MODES = ['major', 'minor']


def next_id(current, mode='minor'):
    if mode == 'none' or mode is None:
        return str(current)

    if mode not in MODES:
        raise UnknownReleaseMode(mode)

    current = parse_id(current)
    if mode == 'minor':
        current['minor'] += 1

    if mode == 'major':
        current['minor'] = 0
        current['major'] += 1

    return format_id(current)


def parse_id(release_id):
    try:
        parts = map(int, release_id.split('.'))
    except:
        raise BadlyFormattedRelease("Can't process release id: %s" %
                                    release_id)

    if len(parts) != 2:
        raise BadlyFormattedRelease("Release id must have major an minor")

    return dict(zip(MODES, parts))


def format_id(current):
    return '.'.join((str(current['major']), str(current['minor'])))
