import sys
import csv

from pymotifs.reports import nr as _nr


def write(headers, rows):
    writer = csv.DictWriter(sys.stdout, headers, delimiter="\t")
    writer.writerow(dict(zip(headers, headers)))
    for row in rows:
        writer.writerow(row)


def nr(**kwargs):
    data = _nr.report(kwargs.pop('version'), kwargs.pop('resolution'), **kwargs)
    write(_nr.HEADERS, data)
