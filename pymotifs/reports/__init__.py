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


def write(headers, rows):
    writer = csv.DictWriter(sys.stdout, headers, delimiter="\t")
    writer.writerow(dict(zip(headers, headers)))
    for row in rows:
        writer.writerow(row)


def nr_groups(**kwargs):
    session = setup(kwargs['engine'])
    version = kwargs.pop('version')
    resolution = kwargs.pop('resolution')
    data = _nr.groups(session, version, resolution, **kwargs)
    write(data.headers, data.rows)


def nr_pairs(**kwargs):
    session = setup(kwargs['engine'])
    version = kwargs.pop('version')
    resolution = kwargs.pop('resolution')
    data = _nr.pairs(session, version, resolution, **kwargs)
    write(data.headers, data.rows)


def species(**kwargs):
    session = setup(kwargs['engine'])
    data = _species.report(session, **kwargs)
    write(_species.HEADERS, data)
