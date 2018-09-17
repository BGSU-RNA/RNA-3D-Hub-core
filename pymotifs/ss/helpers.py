# read a postscript file, extract nucleotide sequence and x,y coordinates, with scale and translation applied
# returns a list of dictionaries, each of which has x_coordinate, y_coordinate, and unit (base)
import re
import hashlib

from pymotifs import core


def md5(filename):
    with open(filename, 'rb') as raw:
        return hashlib.md5(raw.read()).hexdigest()


class Base(object):
    @classmethod
    def matches(cls, line):
        return cls.pattern.match(line)


class Nucleotide(Base):
    pattern = \
        re.compile('^\((?P<sequence>\w)\) (?P<x>-?\d+\.\d+) ' +
                   '(?P<y>-?\d+\.\d+) lwstring$')

    def __init__(self, x=0, y=0, sequence=''):
        self.sequence = sequence.upper()
        if self.sequence == 'Y':
            self.sequence = 'U'

        self.x = float(x)
        self.y = float(y)

    def shift(self, nt):
        self.x = self.x - nt.x
        self.y = nt.y - self.y

    def update(self, trans):
        self.x = (self.x * trans.s_x + trans.t_x)
        self.y = (self.y * trans.s_y + trans.t_y)

    def to_point(self):
        point2pixel = lambda number: 96.0/72.0 * number
        return {
            'x_coordinate': round(point2pixel(self.x), 3),
            'y_coordinate': round(point2pixel(self.y), 3),
            'unit': self.sequence
        }

    def __repr__(self):
        return '<Nucleotide x: %s y: %s sequence: %s>' \
            % (self.x, self.y, self.sequence)


class Transform(Base):
    def __init__(self, **kwargs):
        self.s_x = float(kwargs.get('s_x') or 1.0)
        self.s_y = float(kwargs.get('s_y') or 1.0)
        self.t_x = float(kwargs.get('t_x') or 0.0)
        self.t_y = float(kwargs.get('t_y') or 0.0)

    def update(self, trans):
        self.s_x = self.s_x * trans.s_x
        self.s_y = self.s_y * trans.s_y
        self.t_x = self.t_x + self.s_x * trans.t_x
        self.t_y = self.t_y + self.s_y * trans.t_y

    def __repr__(self):
        return '<Transform t_x: %s t_y: %s s_x: %s s_y: %s>' \
            % (self.t_x, self.t_y, self.s_x, self.s_y)


class Scale(Transform):
    pattern = re.compile('^(?P<s_x>\d+.\d+) (?P<s_y>\d+.\d+) scale$')


class Translate(Transform):
    pattern = re.compile('^(?P<t_x>-?\d+\.\d+) (?P<t_y>-?\d+\.\d+) translate$')


def type_match(line):
    types = [Nucleotide, Scale, Translate]
    for builder in types:
        match = builder.matches(line)
        if match:
            return (builder, builder(**match.groupdict()))
    return (None, None)


class Parser(core.Base):
    def parse(self, raw):
        data = []
        transform = Transform()
        shift = Nucleotide(x=100)

        for line in raw:
            line = re.sub('%.*$', '', line.rstrip())
            match, info = type_match(line)
            if match is None and len(data) < 10:
                data = []
            elif match == Nucleotide:
                info.update(transform)
                if info.x < shift.x:
                    shift.x = info.x
                if info.y > shift.y:
                    shift.y = info.y
                data.append(info)
            elif match and issubclass(match, Transform):
                transform.update(info)

        points = []
        for index, entry in enumerate(data):
            entry.shift(shift)
            point = entry.to_point()
            point['index'] = index
            points.append(point)

        return points

    def __call__(self, filename):
        with open(filename, 'rb') as raw:
            parsed = self.parse(raw)

        return {
            'sequence': ''.join(nt['unit'] for nt in parsed),
            'nts': parsed
        }
