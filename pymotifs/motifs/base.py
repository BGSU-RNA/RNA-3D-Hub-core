import os
import csv

from pymotifs import core
from pymotifs import models as mod
from pymotifs.utils.releases import Release


class CsvLoader(core.MassLoader):
    """ Name of the file to load """
    name = None

    """ Name of the table to write to """
    table = None

    """ Headers for the csv file to parse """
    headers = None

    def __init__(self, *args, **kwargs):
        super(CsvLoader, self).__init__(*args, **kwargs)
        if not self.name:
            raise core.InvalidState("Must specify name")
        if not self.table:
            raise core.InvalidState("Must specify the table name")
        if not self.headers:
            raise core.InvalidState("Must specify file headers")

    def release(self):
        return Release(self.config, self.session.maker).current('motif')

    def base(self, loop_type):
        release = self.release()
        with self.session() as session:
            return session.query(mod.MlReleases).\
                filter_by(id=release, type=loop_type).\
                one().description

    def filename(self, loop_type):
        filename = os.path.join(self.config['locations']['releases_dir'],
                                self.base(loop_type), self.name)
        if not os.path.exists(filename):
            raise core.InvalidState("Cannot find motif file %s" % filename)
        return filename

    def parse(self, handle, release_id):
        reader = csv.DictReader(handle, fieldnames=self.headers,
                                delimiter=',', quotechar='"')
        data = []
        for row in reader:
            row['release_id'] = release_id
            data.append(row)
        return data

    def data(self, *args, **kwargs):
        filename = self.filename()
        release_id = self.release()
        with open(filename, 'rb') as raw:
            parsed = self.parse(raw, release_id)

        return [self.table(**entry) for entry in parsed]
