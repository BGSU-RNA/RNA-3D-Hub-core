#!/usr/bin/env python

import os
import sys

pymotifs = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.append(pymotifs)

from pymotifs.cli.commands import cli


if __name__ == '__main__':
    cli()
