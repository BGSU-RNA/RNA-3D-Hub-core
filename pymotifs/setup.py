import random
import logging
from copy import deepcopy

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from pymotifs import models as mod
from pymotifs import utils as util
from pymotifs import config as conf


class Runnable(object):
    def __init__(self, opts):
        self.logger = logging.getLogger(__name__)
        self.options = deepcopy(opts)
        self.stage = self.options.pop('name')
        self._ids = self.options.pop('ids')
        self.config = conf.load(self.options.pop('config'))
        self._connection = None
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
        from pymotifs.dispatcher import Dispatcher
        return Dispatcher(self.stage, self.config, self.connection,
                          **self.options)

    @property
    def pdbs(self):
        if self._ids:
            return self._ids

        helper = util.pdb.RnaPdbsHelper()
        if self.options.get('all', False):
            self.logger.debug("Getting all PDBs")
            return helper()

        if self.options.pop('known', None):
            self.logger.debug("Using known pdbs only")
            return list(util.known(self.config, pdb=False))

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


def loggers(options):
    opts = deepcopy(options)
    log_args = {
        'level': getattr(logging, opts.pop('log_level', 'info').upper()),
        'filemode': opts.pop('log_mode'),
        'format': '%(levelname)s:%(name)s:%(asctime)s:%(message)s',
    }

    filename = opts.pop('log_file', None)
    if filename:
        log_args['filename'] = filename

    logging.basicConfig(**log_args)
    # logging.captureWarnings(True)
    base = logging.getLogger()
    pool_logger = logging.getLogger('sqlalchemy.pool')
    pool_logger.setLevel(logging.ERROR)
    for handler in base.handlers:
        pool_logger.addHandler(handler)

    return opts


def seed(options):
    if 'seed' in options:
        random.seed(options.pop('seed'))
    return options


def setup(given):
    opts = loggers(given)
    opts = seed(opts)
    return Runnable(opts)
