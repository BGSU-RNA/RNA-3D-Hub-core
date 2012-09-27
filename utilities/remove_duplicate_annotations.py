"""

A helper script used when debugging motif annotation inheritance.

"""

import sys
import os.path

# add parent directory to the path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), os.path.pardir)))

import models


session = models.session
MotifAnnotation = models.MotifAnnotation

for annotation in session.query(MotifAnnotation).all():
    if annotation.annotation is not None:
        new_annotation = ', '.join(set(annotation.annotation.split(', ')))
        if new_annotation != annotation.annotation:
            print '%s vs %s' % (annotation.annotation, new_annotation)
            annotation.annotation = new_annotation

    if annotation.author is not None:
        new_author = ', '.join(set(annotation.author.split(', ')))
        if new_author != annotation.author:
            print '%s vs %s' % (annotation.author, new_author)
            annotation.author = new_author

    if annotation.common_name is not None:
        new_common_name = ', '.join(set(annotation.common_name.split(', ')))
        if new_common_name != annotation.common_name:
            print '%s vs %s' % (annotation.common_name, new_common_name)
            annotation.common_name = new_common_name

    session.merge(annotation)

session.commit()
