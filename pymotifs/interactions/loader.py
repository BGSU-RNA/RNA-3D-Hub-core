"""Run all interaction level stages. This runs:

interactions.pairwise
    Annotate all pairwise interactions
interactions.flanking
    Annotate all flanking interactions
interactions.summary
    Summarize the number of interactions for each unit.
"""

import pymotifs.core as core

from pymotifs.interactions.annotate_python import Loader as AnnotationLoader


class Loader(core.StageContainer):
    stages = set([AnnotationLoader])
