import re
import logging

from sqlalchemy import Table
from sqlalchemy import Column
from sqlalchemy import String
from sqlalchemy import Integer
from sqlalchemy import MetaData

from sqlalchemy.ext.compiler import compiles
from sqlalchemy.sql.expression import Delete

from sqlalchemy.ext.declarative import declarative_base


logger = logging.getLogger(__name__)

metadata = MetaData()
Base = declarative_base(metadata=metadata)


def camelize_classname(tablename):
    return str(tablename[0].upper() +
               re.sub(r'_(\w)', lambda m: m.group(1).upper(), tablename[1:]))


def define_missing_views(metadata):

    # We have to define what columns are PK for views
    Table('exp_seq_pdb', metadata,
          Column('exp_seq_id', String, primary_key=True),
          Column('pdb_id', String, primary_key=True),
          Column('chain_id', Integer, primary_key=True),
          extend_existing=True)

    Table('correspondence_pdbs', metadata,
          Column('correspondence_id', String, primary_key=True),
          Column('chain_id_1', Integer, primary_key=True),
          Column('chain_id_2', Integer, primary_key=True),
          extend_existing=True)


def reflect(engine):
    """Reflect all tables/views from the database into python. This cannot
    reflect any table/view that does not have a primary key as this is a
    limitation of sqlalchemy.

    Modified from
    https://charleslavery.com/notes/sqlalchemy-reflect-tables-to-declarative.html
    """

    metadata.bind = engine
    metadata.reflect(views=True)
    define_missing_views(metadata)

    glo = globals()
    for name, obj in metadata.tables.items():
        classname = camelize_classname(name)

        try:
            glo[classname] = type(classname, (Base,), {'__table__': obj})
        except:
            logger.warning("Could not reflect table %s", name)


@compiles(Delete, 'mysql')
def compile_delete_with_joins(element, compiler, **kwargs):
    """This is added to allow sqlalchemly to compile delete's with join
    statements. This is covered in:

    https://bitbucket.org/zzzeek/sqlalchemy/issues/959/support-mysql-delete-from-join

    when that gets fixed we can remove this code.The idea behind this is to use
    some mysql specific syntax to allow for delete's with join statements. The
    normal compilation does not allow for join statements in delete. This
    should.
    """

    if element.table._is_join:
        name = element.table.left.name
    else:
        name = element.table.name

    table = compiler.process(element.table, asfrom=True)

    text = 'DELETE %s FROM %s'
    text = text % (name, table)

    if element._whereclause is not None:
        text += ' WHERE ' + compiler.process(element._whereclause)

    return text
