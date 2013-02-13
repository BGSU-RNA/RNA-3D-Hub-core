"""


"""
import re

from MotifAtlasBaseClass import MotifAtlasBaseClass
from models import session, MotifAnnotation


for loop in session.query(MotifAnnotation).all():

    if loop.common_name is not None and loop.common_name != '':
        loop.common_name = re.sub('\s?\[.+?\]\s?', '', loop.common_name)
        loop.common_name = ' | '.join(set(loop.common_name.split('| ')))

    if loop.annotation is not None and loop.annotation != '':
        loop.annotation = re.sub('\s?\[.+?\]\s?', '', loop.annotation)
        loop.annotation = ' | '.join(set(loop.annotation.split('| ')))



session.commit()