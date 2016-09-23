import sys
import csv

from sqlalchemy.orm import sessionmaker

from pymotifs.core import Session
from pymotifs import models as mod

from pymotifs.reports import nr as _nr
from pymotifs.reports import species as _species


def setup(engine):
    mod.reflect(engine)
    return Session(sessionmaker(engine))


def write(headers, rows, hide_headers=True, delimiter='\t', **kwargs):
    writer = csv.DictWriter(sys.stdout, headers, delimiter=str(delimiter))
    if not hide_headers:
        writer.writerow(dict(zip(headers, headers)))
    for row in rows:
        writer.writerow(row)


def nr_groups(**kwargs):
    config = kwargs.pop('config')
    version = kwargs.pop('version')
    resolution = kwargs.pop('resolution')
    stage = _nr.Groups(config, setup(kwargs['engine']))
    stage((version, resolution), **kwargs)


def nr_pairs(**kwargs):
    config = kwargs.pop('config')
    version = kwargs.pop('version')
    resolution = kwargs.pop('resolution')
    stage = _nr.Pairs(config, setup(kwargs['engine']))
    stage((version, resolution), **kwargs)


def species(**kwargs):
    session = setup(kwargs['engine'])
    data = _species.report(session, **kwargs)
    write(_species.HEADERS, data, **kwargs)
