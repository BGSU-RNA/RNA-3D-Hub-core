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


class KeyValueType(click.ParamType):
    name = 'key_value'

    def convert(self, value, param, ctx):
        if value is None:
            return None
        if re.match('\w+=\w+', value):
            return value
        self.fail("Invalid key value pair")

    def finalize_dict(self, value):
        data = {}
        for entry in value:
            key, value = entry.split('=')
            data[key] = value
        return data


PDB = PdbParamType()
DATE = DateParamType()
KEY_VALUE = KeyValueType()
FILE = click.Path(exists=True, dir_okay=False, resolve_path=True)
