import random
import logging
import collections as coll

import click

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from pymotifs.cli import setup
from pymotifs.cli import introspect
from pymotifs.cli.params import FILE
from pymotifs.cli.params import DATE
from pymotifs.cli.params import PDB
from pymotifs.constants import BOOTSTRAPPING_PDBS
from pymotifs.email import Emailer

from pymotifs import models as mod
from pymotifs import config as conf
from pymotifs.version import __VERSION__
from pymotifs.dispatcher import Dispatcher


# This constant is included and then changed during testing so we do not
# actually execute the stages unless we need to. This allows us to better test
# if we are loading all the options and configuration correctly.
SHOULD_RUN = True


Result = coll.namedtuple('Result', ['result', 'success', 'message'])


class Runnable(object):
    """This object works to define the what will be run. This is mostly used to
    make it easier to inspect and test out commands. I've found a couple issues
    where things aren't be used correctly so I'm restructing things to create
    more tests.
    """

    def __init__(self, name, ids, config, engine, options):
        self.config = config
        self.target = name
        self._ids = ids
        self.engine = engine
        self._options = options

    @property
    def skipped_stages(self):
        return setup.expand_stage_pattern(self.target,
                                          self._options.get('skip_stage'))

    @property
    def recalculate_stages(self):
        return setup.expand_stage_pattern(self.target,
                                          self._options.get('recalculate'))

    @property
    def session(self):
        return sessionmaker(self.engine)

    @property
    def dispatcher(self):
        return Dispatcher(self.target, self.config, self.session,
                          skip_dependencies=self.skip_dependencies,
                          exclude=self.skipped_stages)

    @property
    def skip_dependencies(self):
        return self.options.get('skip_dependencies', False)

    @property
    def options(self):
        opts = dict(self._options)
        opts['recalculate'] = self.recalculate_stages
        opts['dry_run'] = SHOULD_RUN and opts['dry_run']
        return opts

    def ids(self):
        return self._ids or setup.pdbs(self.config, self.options)

    def __call__(self):
        if 'seed' in self.options:
            random.seed(self.options['seed'])

        result = None
        success = None
        try:
            if SHOULD_RUN:
                result = self.dispatcher(self.ids(), **self.options)
            success = True
        except Exception as error:
            success = False
            logging.exception(error)

        msg = None
        if self.options['email']:
            emailer = Emailer(self.config, self.engine)
            emailer = emailer.success if success else emailer.fail
            msg = emailer(self.target,
                          log_file=self.options['log_file'],
                          dry_run=self.options['dry_run'])

        return Result(result=result, success=success, message=msg)


def run(ctx, name, ids, config=None, engine=None, **kwargs):
    mod.reflect(engine)
    runnable = Runnable(name, ids, config, engine, kwargs)
    result = runnable()
    if not result.success:
        click.secho("Pipeline failed", fg='red', err=True)
        ctx.exit(1)


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
        'engine': create_engine(config['db']['uri'])
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
              type=DATE, help='Get PDBs from after the date')
@click.option('--before-date', default=None,
              type=DATE, help='Get PDBs from before the date')
@click.option('--exclude', multiple=True, type=PDB,
              help='Excluded PDB(s)')
@click.option('--ignore-time', is_flag=True, help='Ignore time for rerunning')
@click.option('--use-pdb', is_flag=True, help='Use PDB files instead of cif')
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


@cli.command(short_help='populate a database')
@click.pass_context
def bootstrap(ctx, **kwargs):
    """Poupulate a testing database.

    To test this pipeline working copy of the database is needed. This command
    will populate a database with some default data for testing. This will not
    import any distance data.
    """

    kwargs['skip_stage'] = ['units.distances', 'export', 'pdb.obsolete']
    kwargs['seed'] = 1
    kwargs['config'] = 'conf/bootstrap.json'
    kwargs.update(ctx.parent.objs)
    run(ctx, 'update', BOOTSTRAPPING_PDBS, **kwargs)


@cli.command(short_help="List stages")
@click.pass_context
def list(ctx, **kwargs):
    """List all known stages.
    """
    mod.reflect(ctx.parent.objs['engine'])

    formatter = click.HelpFormatter(width=90)
    formatter.write_dl((s[0], s[1]) for s in introspect.stages())
    click.echo(formatter.getvalue(), nl=False)


@cli.command(short_help="Get information about a stage")
@click.argument('name', type=str)
@click.pass_context
def about(ctx, name=None, **kwargs):
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
    """Commands dealing with importing 2D diagrams
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

    This will determine the expermental sequence of the given chain and then
    align it to the secondary structure in the given file. The 2D will be
    imported as needed. Doing one such import can allow the inference of
    positions across many structures since experimental sequences are often
    reused.
    """
    kwargs.update(ctx.parent.objs)
    run(ctx, 'ss.position_mapping', [pdb], **kwargs)
