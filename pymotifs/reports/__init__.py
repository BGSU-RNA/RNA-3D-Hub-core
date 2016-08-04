import sys
import csv

from sqlalchemy.orm import sessionmaker

from pymotifs.core import Session
from pymotifs import models as mod

from pymotifs.reports import nr as _nr
from pymotifs.reports.setup import setup
from pymotifs.reports import species as _species


def write(headers, rows):
    writer = csv.DictWriter(sys.stdout, headers, delimiter="\t")
    writer.writerow(dict(zip(headers, headers)))
    for row in rows:
        writer.writerow(row)

def nr(**kwargs):
    session = setup(kwargs['engine'])
    version =  kwargs.pop('version')
    resolution = kwargs.pop('resolution')
    data = _nr.report(session, version, resolution, **kwargs)
    write(_nr.HEADERS, data)


def species(**kwargs):
    session = setup(kwargs['engine'])
    data = _species.report(session, **kwargs)
    write(_species.HEADERS, data)
