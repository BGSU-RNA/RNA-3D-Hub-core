import os
import tempfile
import pickle


def filename(group):
    """Get a filename for temprorary data. This will create one temprorary file
    per group that may or may not be accesible between runs.

    :group: The name of the group to use
    :returns: The temporary filename
    """

    return os.path.join(tempfile.gettempdir() + '.temp')


def cleanup(group):
    """Remove the temp file for the group.

    :group: The name of the group.
    """

    name = filename(group)
    if os.path.exists(name):
        os.remove(name)
    return True


def store(group, data):
    with open(filename(group), 'wb') as raw:
        pickle.dump(data, raw)


def load(group):
    name = filename(group)
    if not os.path.exists(name):
        return None

    with open(name, 'rb') as raw:
        return pickle.load(raw)
