import re
from datetime import datetime

import click


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


PDB = PdbParamType()
DATE = DateParamType()
FILE = click.Path(exists=True, dir_okay=False, resolve_path=True)
