import os
import re
import inspect
from datetime import datetime

import click

from pymotifs.core import Stage


class DateParamType(click.ParamType):
    name = 'date'

    def convert(self, value, param, ctx):
        try:
            return datetime.strptime(value, '%Y-%m-%d').date()
        except:
            self.fail("Invalid date, must by y-m-d")


class PdbParamType(click.ParamType):
    name = 'pdb'

    def convert(self, value, param, ctx):
        if re.match('^[1-9][0-9a-z]{3}$', value, re.IGNORECASE):
            return value
        self.fail("Invalid PDB id")


class StageCommand(click.Command):
    def __init__(self, name, **kwargs):
        kwargs['callback'] = self.run_stage
        super(name, **kwargs)
        import_name = 'pymotifs.' + name
        from_list = import_name.split('.')
        module = __import__(import_name, fromlist=from_list)
        pairs = inspect.getmembers(module, self.is_loader(import_name))
        if not pairs:
            raise ImportError("Could not find a class to use from: %s" % name)

        self.stage = pairs[0][1]

    def run_stage(self, **kwargs):
        pass


class StageGroup(click.MultiCommand):

    def get_loaders(self, base):
        loaders = []
        ignore = set(['utils', 'cli', 'core', '__init__.py', 'dispatcher.py'])
        for path in os.listdir(base):
            if path.endswith('.py') and path not in ignore:
                name = path[0:-3]
                module = __import__(base + '.' + path[0:-3])
                if bool(inspect.getmembers(module, self.__stage_checker__)):
                    loaders.append(name)
        return loaders

    def __stage_checker__(self, obj):
        return inspect.isclass(obj) and issubclass(obj, Stage) and obj != Stage

    def list_commands(self, ctx):
        return self.get_loaders('pymotifs')

    def get_command(self, ctx, name):
        command = StageCommand(name)
        dec = click.argument('ids', type=PDB)
        return dec(command)


PDB = PdbParamType()
DATE = DateParamType()
FILE = click.Path(exists=True, dir_okay=False, resolve_path=True)
