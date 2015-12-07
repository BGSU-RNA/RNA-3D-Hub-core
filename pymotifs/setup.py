import logging

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from pymotifs.utils import known
from pymotifs import models as mod
from pymotifs import config as conf
from pymotifs.dispatcher import Dispatcher
from pymotifs.utils.pdb import RnaPdbsHelper


def log(options):
    log_args = {
        'level': getattr(logging, options['log_level'].upper()),
        'filemode': options['log_mode'],
        'format': '%(levelname)s:%(name)s:%(asctime)s:%(message)s',
    }

    if 'log_file' in options:
        log_args['filename'] = options['log_file']

    logging.basicConfig(**log_args)
    # logging.captureWarnings(True)
    base = logging.getLogger()
    pool_logger = logging.getLogger('sqlalchemy.pool')
    pool_logger.setLevel(logging.ERROR)
    for handler in base.handlers:
        pool_logger.addHandler(handler)


def connection(config):
    engine = create_engine(config['db']['uri'])
    mod.reflect(engine)
    return sessionmaker(bind=engine)


def config(options):
    return conf.load(options['config'])


class Runnable(object):
    def __init__(self, opts):
        self.logger = logging.getLogger(__name__)
        self.options = opts
        self.stage = self.options.pop('name')
        self._ids = self.options.pop('ids')
        self.config = self.options.pop('config')
        self._connection = self.options.pop('connection', None)
        self.__recalculate__()

    @property
    def connection(self):
        if not self._connection:
            engine = create_engine(self.config['db']['uri'])
            mod.reflect(engine)
            self._connection = sessionmaker(bind=engine)
        return self._connection

    @property
    def dispatcher(self):
        return Dispatcher(self.stage, self.config, self.connection,
                          **self.options)

    @property
    def pdbs(self):
        if self._ids:
            return self._ids

        helper = RnaPdbsHelper()
        if self.options.get('all', False):
            self.logger.debug("Getting all PDBs")
            return helper()

        if self.options.pop('known', None):
            self.logger.debug("Using known pdbs only")
            return list(known(self.config, pdb=False))

        kwargs = {}
        if 'before' in self.options or 'after' in self.options:
            kwargs['dates'] = (self.options.get('after', None),
                               self.options.get('before', None))
            self.logger.debug("Geting PDBs within dates %s, %s",
                              *kwargs['dates'])
            return helper(**kwargs)

        return []

    def __recalculate__(self):
        updated = []
        for recalc in self.options.get('recalculate', []):
            if recalc == '.':
                updated.append('pymotifs.' + self.stage)
            elif recalc == '*':
                updated = True
                break
            else:
                updated.append('pymotifs.' + recalc)

        if updated:
            self.options['recalculate'] = updated
