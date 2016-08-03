from sqlalchemy.orm import sessionmaker

from pymotifs.core import Session
from pymotifs import models as mod


def setup(engine):
    mod.reflect(engine)
    return Session(sessionmaker(engine))
