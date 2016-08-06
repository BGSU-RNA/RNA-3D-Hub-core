def units(column, **kwargs):
    mod.reflect(kwargs['engine'])
    table_name, name = column.split('.')
    table_name = mod.camelize_classname(table_name)
    table = getattr(mod, table_name)
    session = Session(sessionmaker(kwargs['engine']))
    corrector = TableCorrector(kwargs['config'], session)
    corrector(table, name, **kwargs)


def nr_history(version, **kwargs):
    mod.reflect(kwargs['engine'])
    session = Session(sessionmaker(kwargs['engine']))
    # TODO: Load the given release into a cached file
    # Run the nr.parent_counts stage with some PDBs make sure to use
    # recalculate=True
