import inspect
import functools as ft

from pymotifs.core import Stage

BASE = 'pymotifs'


def get_stages(name, obj):
    return inspect.isclass(obj) and \
        issubclass(obj, Stage) and \
        obj.__module__ == name and \
        obj != Stage


def imports(seen, obj):
    return inspect.ismodule(obj) and \
        BASE in obj.__name__ and \
        obj not in seen


def has_stages(module):
    isstage = ft.partial(get_stages, module.__name__)
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
    parts = name.split('.')
    fromlist = [BASE]
    fromlist.extend(parts[:-1])
    return __import__(BASE + '.' + name, fromlist=fromlist)


def get_stage_info(name):
    module = get_stage(name)
    return stage_info(module)


def get_loader(name):
    module = get_stage(name)
    isstage = ft.partial(get_stages, module.__name__)
    loaders = inspect.getmembers(module, isstage)
    return next(loaders)


def traverse(module, seen):
    seen.add(module)

    stages = set()
    if has_stages(module):
        stages.add(stage_info(module))

    isimport = ft.partial(imports, seen)
    for name, includes in inspect.getmembers(module, isimport):
        stages.update(traverse(includes, seen))

    return stages


def stages():
    update = __import__(BASE + '.update')
    return sorted(traverse(update, set()))
