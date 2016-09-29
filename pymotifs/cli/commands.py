"""This contains all command line interface commands that the pipeline script
can run. This contains all the functions that will be run.
"""

import random
import logging

import click

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from pymotifs.cli import setup
from pymotifs.cli import introspect
from pymotifs.cli.params import FILE
from pymotifs.cli.params import DATE
from pymotifs.cli.params import PDB
from pymotifs.constants import BOOTSTRAPPING
from pymotifs.email import Emailer

from pymotifs import correct as _correct
from pymotifs import models as mod
from pymotifs import config as conf
from pymotifs import transfer as _transfer
from pymotifs.version import __VERSION__
from pymotifs.dispatcher import Dispatcher


def run(ctx, name, ids, config=None, engine=None, **kwargs):
    """Actually run the pipeline. This is the main function which will load and
    run the stages requested. It fetch PDBs to run if needed, reflect the
    models, run stages, and send an email as needed.

    Parameters
    ----------
    name : str
        The name of the stage to run.
    ids : list
        The PDB ids to use.
    config : dict
        The configuration dictionary as produced by pymotifs.config.load.
    engine : Engine
    The SQL Alchemy engine connection.
    **kwargs : dict, optional
        The other keyword arguments which will be passed to dispatcher.
    """

    if kwargs.get('seed', None) is not None:
        random.seed(kwargs['seed'])

    mod.reflect(engine)

    try:
        setup.expand_stage_pattern(name, 'recalculate', kwargs)
        setup.expand_stage_pattern(name, 'skip_stage', kwargs)
    except introspect.UnknownStageError as err:
        click.secho("Unknown stage %s" % err.args, err=True, fg='red')
        ctx.exit(1)

    if not ids:
        ids = setup.pdbs(config, kwargs)

    kwargs['exclude'] = kwargs.get('skip_stage')

    error = None
    dispatcher = Dispatcher(name, config, sessionmaker(engine), **kwargs)
    mailer = Emailer(config, engine)
    try:
        dispatcher(ids, **kwargs)
    except Exception as error:
        click.secho("Pipeline failed", fg='red', err=True)
        logging.exception(error)
        ctx.exit(1)
    finally:
        if kwargs['email']:
            mailer(name, ids=ids, error=None, **kwargs)


@click.group(short_help="Interface to update pipeline",
             context_settings={'help_option_names': ['-h', '--help']})
@click.option('--config', default='conf/motifatlas.json', type=FILE,
              help="JOSN configuration file")
@click.option('--log-file', type=click.Path(dir_okay=False, resolve_path=True),
              help="Logging filename")
@click.option('--log-level', default='info',
              type=click.Choice(['debug', 'info', 'warning', 'error']),
              help="Logging level to use")
@click.option('--log-mode', default='a', type=click.Choice(['w', 'a']),
              help='Mode to open the  logging file')
@click.option('--email/--no-email', default=False, help='Send email')
@click.version_option(__VERSION__)
@click.pass_context
def cli(ctx, **options):
    """A tool for running parts of the update pipeline.

    This allows users run either the whole or parts of the update pipeline. It
    offers fine grained control over what gets run and how.
    """
    setup.logs(options)
    ctx.objs = options
    config = conf.load(options['config'])
    ctx.objs.update({
        'config': config,
        'engine': create_engine(config['db']['uri'],
                                pool_size=config['db'].get('pool_size', 20),
                                max_overflow=config['db'].get('max_overflow', 0))
    })


@cli.command('run', short_help="Run stage(s)")
@click.option('--dry-run', is_flag=True, help="Alter nothing while running")
@click.option('--skip-dependencies', is_flag=True, help='Skip stage(s)')
@click.option('--skip-stage', multiple=True, help='Stage to skip')
@click.option('--seed', type=int, help="Set the random seed")
@click.option('--recalculate', multiple=True,
              help="Recalculate data for the given stage(s)")
@click.option('--all', is_flag=True, help="Use all RNA containing PDBS")
@click.option('--known', is_flag=True, help="Use only downloaded files")
@click.option('--after-date', default=None,
              type=DATE, help='Get PDBs from after the date (YYYY-MM-DD)')
@click.option('--before-date', default=None,
              type=DATE, help='Get PDBs from before the date (YYYY-MM-DD)')
@click.option('--exclude', multiple=True, type=PDB,
              help='Excluded PDB(s)')
@click.option('--ignore-time', is_flag=True, help='Ignore time for rerunning')
@click.argument('name')
@click.argument('ids', nargs=-1, type=PDB)
@click.pass_context
def do(ctx, name, ids, **kwargs):
    """Run specific stages in the pipeline.

    This accepts a stage to run and a list of PDB ids to import. It determines
    which if any dependencies must be run and runs them all.
    """
    kwargs.update(ctx.parent.objs)
    run(ctx, name, ids, **kwargs)


@cli.command(short_help='Populate a database')
@click.pass_context
def bootstrap(ctx, **kwargs):
    """Populate a testing database.

    To test this pipeline working copy of the database is needed. This command
    will populate a database with some default data for testing. This will not
    export interactions or loops, nor does it load the obsolete data. This will
    run 4 updates to create several nr and motif data sets to use.
    """

    kwargs.update(ctx.parent.objs)
    kwargs['config'] = conf.load('conf/bootstrap.json')
    kwargs['engine'] = create_engine(kwargs['config']['db']['uri'])
    for index, step in enumerate(BOOTSTRAPPING['steps']):
        logging.info("Running step %i", index)
        args = dict(kwargs)
        args.update(BOOTSTRAPPING['args'])
        args.update(step.get('args', {}))
        run(ctx, step['stage'], step['pdbs'], **args)


@cli.command(short_help="List stages")
@click.pass_context
def list(ctx, **kwargs):
    """This will output each stage name along with the short documentation on
    the stage, if any is available. This will only list stages which are part
    of the 'update' stage.
    """

    mod.reflect(ctx.parent.objs['engine'])
    formatter = click.HelpFormatter(width=90)
    formatter.write_dl((s[0], s[1]) for s in introspect.stages())
    click.echo(formatter.getvalue(), nl=False)


@cli.command(short_help="Get information about a stage")
@click.argument('name', type=str)
@click.pass_context
def about(ctx, name=None, **kwargs):
    """Display help for a stage. This will write out the complete documentation
    for each stage with the given name.
    """

    mod.reflect(ctx.parent.objs['engine'])

    info = introspect.get_stage_info(name)
    if not info:
        click.echo("Unknown stage %s" % name, err=True, fg='red')
        ctx.exit(1)

    formatter = click.HelpFormatter()
    with formatter.section(name=name):
        formatter.write_text(info[2])
    click.echo(formatter.getvalue())


@cli.group(short_help="Run commands for 2d diagrams")
@click.pass_context
def ss(ctx):
    """Commands dealing with importing 2D diagrams.
    """
    ctx.objs = ctx.parent.objs


@ss.command('import', short_help='Import a 2D diagram')
@click.option('--ss-name', default=None, type=str, help='Name for the diagram')
@click.option('--recalculate', multiple=True,
              help="Recalculate data for the given stage(s)")
@click.argument('filename', metavar='F', type=FILE)
@click.pass_context
def ss_import(ctx, **kwargs):
    """Import a 2D diagram.

    This will parse the secondary structure diagram in the given file and
    import it. The file should be a post script file from the gutell lab. It
    parses their style of postcript files and no others.
    """
    kwargs.update(ctx.parent.objs)
    run(ctx, 'ss.positions', None, **kwargs)


@ss.command('align', short_help='Align a 2D to a chain')
@click.option('--ss-name', default=None, type=str, help='Name for the diagram')
@click.option('--recalculate', multiple=True,
              help="Recalculate data for the given stage(s)")
@click.option('--skip-dependencies', is_flag=True,
              help='Skip running any dependencies')
@click.argument('pdb', type=PDB)
@click.argument('chain')
@click.argument('filename', type=FILE)
@click.pass_context
def ss_align(ctx, pdb, **kwargs):
    """Align a chain to a secondary structure.

    This will determine the experimental sequence of the given chain and then
    align it to the secondary structure in the given file. The 2D will be
    imported as needed. Doing one such import can allow the inference of
    positions across many structures since experimental sequences are often
    reused.
    """
    kwargs.update(ctx.parent.objs)
    run(ctx, 'ss.position_mapping', [pdb], **kwargs)

@cli.group(short_help="A group of commands for correcting issues in the db")
@click.pass_context
def correct(ctx):
    """Correct issues in the database.

    During the migration in the format of the database there were a variety of
    issues found and raised. For example, it was found that not all old style
    ids were translated to correct new style. This group of commands exists to
    make correcting these issues easy.
    """
    ctx.objs = ctx.parent.objs


@correct.command('units', short_help='Correct unit ids in some table')
@click.option('--size', default=1000, type=int,
              help='Number of rows to change at once')
@click.option('--translate', default=False, is_flag=True,
              help='If unit ids should be translated to new style as well')
@click.argument('column', type=str)
@click.pass_context
def correct_units(ctx, column, **kwargs):
    """Correct a column of unit ids in the database.

    This will go through a column in the database and attempt to correct the
    unit ids. This is useful because some old style ids were translated in
    correctly. In addition, we can also try to translate old style if needed.
    This will fail if any unit id in the column is unfixable. It can also
    handle a column that is a comma seperated string of unit ids.

    Running translation requires that the __pdb_unit_id_correspondence table
    exists. Running correction requires that the unit_info table is populated
    with all units from all structures that will be corrected.
    """
    kwargs.update(ctx.parent.objs)
    _correct.units(column, **kwargs)


@correct.command('nr', short_help='Correct NR history counts')
@click.argument('version', type=str)
@click.pass_context
def correct_nr_history(ctx, version, **kwargs):
    """Correct NR history counts.

    When we transitioned the logic of building NR sets, we messed up the logic
    for getting the counts of changes correctly. This recomputes the counts for
    a given NR release.
    """
    pass


@cli.group(short_help='Create reports')
@click.option('--hide-headers', is_flag=True, default=False,
              help="Hide headers in output")
@click.option('--delimiter', default='\t', type=str,
              help='Delimeter to separate  columns with')
@click.pass_context
def report(ctx, **kwargs):
    """Commands dealing with creating reports
    """
    ctx.objs = ctx.parent.objs
    ctx.objs.update(kwargs)
    base = logging.getLogger()
    base.setLevel(logging.CRITICAL)


@report.group('nr', short_help='Reports about the NR set')
@click.pass_context
def report_nr(ctx):
    """Create reports about NR groups
    """
    ctx.objs = ctx.parent.objs


@report_nr.command('groups', short_help='Create a report for an nr set')
@click.option('--resolution', default='all', type=str,
              help='The resolution cutoff to use')
@click.argument('version')
@click.pass_context
def report_nr_groups(ctx, **kwargs):
    """Create a report of the NR set.

    This will detail the members of each set and information about them. It
    creates a CSV that is written to stdout. This report is meant to be useful
    for evaluating the grouping of an NR set.
    """
    kwargs.update(ctx.parent.objs)
    ids = (kwargs.pop('version'), kwargs.pop('resolution'))
    run(ctx, 'reports.nr.groups', ids, **kwargs)


@report_nr.command('pairs',
                   short_help='Create a report about pairs in EQ sets')
@click.option('--resolution', default='all', type=str,
              help='The resolution cutoff to use')
@click.argument('version')
@click.pass_context
def report_nr_pairs(ctx, **kwargs):
    """Create a report about pairs in the NR set.

    This will create a report about the alignment and discrepancy values for
    all pairs within all equivleance sets in the requested NR set. Singleton
    groups do not appear because they have no pairs within them.

    This will produce a CSV that is written to stdout.
    """
    kwargs.update(ctx.parent.objs)
    ids = (kwargs.pop('version'), kwargs.pop('resolution'))
    run(ctx, 'reports.nr.pairs', ids, **kwargs)


@report_nr.command('representative',
                   short_help='Create a report about representative selection')
@click.option('--resolution', default='all', type=str,
              help='The resolution cutoff to use')
@click.argument('version')
@click.pass_context
def report_nr_rep(ctx, **kwargs):
    """Create a report about the selection of representative in the NR set.

    This will create a report about the selection of the representative in the
    NR set. It will write out all members of all groups selected and show the
    number of BP, and NT that were considered in selecting the representative.
    """
    kwargs.update(ctx.parent.objs)
    ids = (kwargs.pop('version'), kwargs.pop('resolution'))
    run(ctx, 'reports.nr.bp_nt', ids, **kwargs)


@report.command('species',
                short_help='Create report about species assignments')
@click.pass_context
def report_species(ctx, **kwargs):
    """Create a report about the species assignments.

    This will look for sequences assigned to species that seem unlikely or are
    assigned to more than one species. It attempts to find assignments that
    seem unlikely.
    """
    kwargs.update(ctx.parent.objs)
    run(ctx, 'reports.species', [], **kwargs)
    # reports.species(**kwargs)


@cli.group(short_help='Dump/Import data')
@click.pass_context
def transfer(ctx):
    """Make dumping data from one database to another easier. This can't easily
    be done using something like mysql-dump because several tables use integer
    primary keys which may be different between different versions of the
    database. This will dump data out and translate the ids as needed.
    """
    ctx.objs = ctx.parent.objs


@transfer.group('cc', short_help='Dump/Import chain chain data')
@click.pass_context
def transfer_chain_chain(ctx):
    """A group of commands for dealing export/import of chain-chain comparisons.
    """
    ctx.objs = ctx.parent.objs


@transfer_chain_chain.command('export', short_help='Dump chain chain data')
@click.argument('filename')
@click.pass_context
def transfer_chain_chain_dump(ctx, **kwargs):
    """Export the chain chain data for import later. Chain-chain comparisions
    can take a long time to compute and so sometimes we want to transfer them
    between machines. This dumps the data into the given file for future
    import.
    """
    kwargs.update(ctx.parent.objs)
    _transfer.chain_chain.dump(**kwargs)


@transfer_chain_chain.command('import', short_help='Import chain chain data')
@click.option('--ignore-missing', is_flag=True,
              help="Ignore any chains that are missing from this DB")
@click.argument('filename')
@click.pass_context
def transfer_chain_chain_import(ctx, **kwargs):
    """This imports exported chain chain comparison data. This will attempt to
    import all data in the given dump file, while ignoring data that already
    exists. It will not error check to see if the database and dump file agree
    on an already existing discrepancy. If it cannot find all chains that are
    in the dump file it will fail, unless the ignore-missing flag is given.
    This is useful when importing a dump from a database that is newer than the
    one that is being imported to.
    """
    kwargs.update(ctx.parent.objs)
    _transfer.chain_chain.load(**kwargs)
