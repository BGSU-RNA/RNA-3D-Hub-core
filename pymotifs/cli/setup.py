import logging

from pymotifs.cli import introspect
from pymotifs.utils import known
from pymotifs.utils.pdb import RnaPdbsHelper


def logs(options):
    log_args = {
        'level': getattr(logging, options['log_level'].upper()),
        'filemode': options['log_mode'],
        'format': '%(levelname)s:%(asctime)s:%(name)s:%(message)s',
    }

    # log_args = {
    #     'level': getattr(logging, options['log_level'].upper()),
    #     'filemode': options['log_mode'],
    #     'format': '%(levelname)s:%(name)s:%(asctime)s:%(message)s',
    # }

    if 'log_file' in options:
        log_args['filename'] = options['log_file']

    logging.basicConfig(**log_args)
    # logging.captureWarnings(True)
    base = logging.getLogger()
    pool_logger = logging.getLogger('sqlalchemy.pool')
    pool_logger.setLevel(logging.ERROR)
    for handler in base.handlers:
        pool_logger.addHandler(handler)


def pdbs(config, options):
    logger = logging.getLogger(__name__)
    helper = RnaPdbsHelper()

    # this is the place to focus on RNA or DNA
    print(options)

    # while developing DNA equivalence classes, run RNA and DNA separately
    nr_molecule_parent_current = config.get('nr_molecule_parent_current','')
    if nr_molecule_parent_current and 'dna' in nr_molecule_parent_current.lower():
        molecule_types = ['DNA',"NA-hybrid"]
    else:
        molecule_types = ['RNA',"NA-hybrid"]
    logger.info("Molecule types: %s" % molecule_types)

    # Note: before-date in the command line becomes before_date in options
    # Note: specifying a before-date or after-date overrides the -all option
    if 'before_date' in options or 'after_date' in options:
        dates = (options.get('after_date', None),
                 options.get('before_date', None))
        logger.info("Getting PDBs within dates %s, %s" % dates)
        return helper(dates=dates,molecule_types=molecule_types)

    if options.get('all', False):
        logger.info("Getting all PDBs")
        return helper(molecule_types=molecule_types)

    if options.pop('known', None):
        logger.info("Using known pdbs only")
        return list(known(config, pdb=False))

    return []


def expand_stage_pattern(stage, key, options):
    updated = []
    for name in options.get(key, []):
        if name == '.':
            name = stage
        elif name == '*':
            updated = True
            break

        if not introspect.is_stage(name):
            raise introspect.UnknownStageError(name)

        updated.append(name)

    if updated:
        options[key] = updated
