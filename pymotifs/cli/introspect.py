"""A module to deal with getting the loaders and things from names. We need to
get a loader to run given it's name as well as getting some documentation. This
is a centralized place for all such activities.
"""

import inspect
import functools as ft

from pymotifs.core import Stage

"""Base name of the module to look up in"""
BASE = 'pymotifs'


class UnknownStageError(Exception):
    """Raised when we have been given a name of stage that does not exist and
    are asked to use it for something.
    """
    pass


def is_stage_loader(module, obj):
    """Useful as a filter for inspect.getmembers for finding all objects which
    are loaders for the current module.
    """
    return inspect.isclass(obj) and \
        issubclass(obj, Stage) and \
        obj.__module__ == module.__name__ and \
        obj != Stage


def is_import(seen, obj):
    """A filter for inspect.getmembers that will pull out modules in the BASE
    module that are not in seen.
    """
    return inspect.ismodule(obj) and \
        BASE in obj.__name__ and \
        obj not in seen


def has_stage_loader(module):
    """Check if the module has a loader.

    :param module: The module to test.
    """
    isstage = ft.partial(is_stage_loader, module)
    return bool(inspect.getmembers(module, isstage))


def stage_info(module):
    name = module.__name__.replace(BASE + '.', '')
    full_doc = inspect.getdoc(module) or ''
    short_doc = ''
    if full_doc:
        end = 30
        if '.' in full_doc:
            end = full_doc.index('.')
        short_doc = full_doc[0:end]
    return (name, short_doc, full_doc)


def get_stage(name):
    """Get the stage with the given name.
    """
    parts = name.split('.')
    fromlist = [BASE]
    fromlist.extend(parts[:-1])
    try:
        return __import__(BASE + '.' + name, fromlist=fromlist)
    except ImportError:
        raise UnknownStageError(name)


def is_stage(name):
    """Check if there is a stage with the given name.
    """
    module = get_stage(name)
    return has_stage_loader(module)


def get_stage_info(name):
    """Loads the information about a stage given it's name.

    :param str name: The name of the stage.
    """
    module = get_stage(name)
    return stage_info(module)


def get_loader(name):
    """Given a name get the loader class for that stage.

    :param str name: Name of the stage to get a loader for.
    """
    module = get_stage(name)
    isstage = ft.partial(is_stage_loader, module)
    loaders = inspect.getmembers(module, isstage)
    if not loaders:
        raise UnknownStageError(name)
    return loaders[0][1]


def traverse(module, seen):
    """Go over all modules imported from a given module and all they import to
    look up all stages. This is how we list all known stages in the pipeline.

    :param module: The module to start with.
    :param seen: A set of already seen modules.
    """
    seen.add(module)

    stages = set()
    if has_stage_loader(module):
        stages.add(stage_info(module))

    isimport = ft.partial(is_import, seen)
    for name, includes in inspect.getmembers(module, isimport):
        stages.update(traverse(includes, seen))

    return stages


def stages():
    """Get all stages part of the main update pipeline.

    :returns: A sorted list of all stages in the update pipeline.
    """
    update = __import__(BASE + '.update')
    stages = sorted(traverse(update, set()))
    return [stage for stage in stages if 'core' not in stage[0]]
