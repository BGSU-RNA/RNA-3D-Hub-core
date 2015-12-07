import random

import click

from pymotifs import setup
from pymotifs.cli.params import FILE
from pymotifs.cli.params import DATE
from pymotifs.cli.params import PDB


@click.group(short_help="Interface to update pipeline",
             context_settings={'help_option_names': ['-h', '--help']})
@click.option('--config', default='conf/motifatlas.json', type=FILE,
              help="JOSN configuration file")
@click.option('--log-file', type=click.Path(dir_okay=False, resolve_path=True),
              help="Logging filename")
@click.option('--log-level', default='debug',
              type=click.Choice(['debug', 'info', 'warning', 'error']),
              help="Logging level to use")
@click.option('--log-mode', default='a', type=click.Choice(['w', 'a']),
              help='Mode to open the  logging file')
@click.option('--email/--no-email', default=False, help='Send email')
@click.version_option('2.0.0')
@click.pass_context
def cli(ctx, **options):
    """A tool for running parts of the update pipeline.
    """
    setup.log(options)
    ctx.objs = options
    config = setup.config(options)
    ctx.objs.update({
        'config': config,
        'connection': setup.connection(config)
    })
    return ctx


@cli.command(short_help="Run stage(s)")
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
@click.argument('name')
@click.argument('ids', nargs=-1, type=PDB)
@click.pass_context
def run(ctx, **kwargs):
    """Run specific stages in the pipeline.
    """
    if kwargs['seed'] is not None:
        random.seed(kwargs['seed'])

    kwargs.update(ctx.parent.objs)
    runnable = setup.Runnable(kwargs)
    runnable.dispatcher(runnable.pdbs, **runnable.options)


@cli.command(short_help='Populate a database')
# bootstraper.set_defaults(config='conf/bootstrap.json')
@click.pass_context
def bootstrap(ctx, **kwargs):
    kwargs['exclude'] = kwargs.get('exclude', [])
    kwargs['exclude'].append('units.distances')
    kwargs['seed'] = 1
    kwargs['name'] = 'update'
    run(ctx, **kwargs)


@cli.group(short_help="Run commands for 2d diagrams")
def ss():
    """Commands dealing with importing 2D diagrams
    """
    pass


@ss.command('import', short_help='Import a 2D diagram')
@click.option('--ss-name', default=None, type=str, help='Name for the diagram')
@click.option('--recalculate', multiple=True,
              help="Recalculate data for the given stage(s)")
@click.argument('filename', metavar='F', type=FILE)
def ss_import(**kwargs):
    kwargs['name'] = ['ss.positions']
    kwargs['ids'] = []
    run(**kwargs)


@ss.command('align', short_help='Align a 2D to a chain')
@click.option('--ss-name', default=None, type=str, help='Name for the diagram')
@click.option('--recalculate', multiple=True,
              help="Recalculate data for the given stage(s)")
@click.option('--skip-dependencies', is_flag=True,
              help='Skip running any dependencies')
@click.argument('pdb', type=PDB)
@click.argument('chain')
@click.argument('filename', type=FILE)
def ss_align(**kwargs):
    kwargs['ids'] = [kwargs.pop('pdb')]
    kwargs['name'] = 'ss.position_mapping'
    run(**kwargs)
