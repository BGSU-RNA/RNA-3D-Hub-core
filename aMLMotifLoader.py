"""

"""

import aMLCollection as motifs
from aMLSqlAlchemyClasses import *


    
A = Uploader(collections=motifs.MotifCollectionMerger(
                    motifs.MotifCollection(file='input/motifs.csv'),
                    motifs.MotifCollection(collection=LoopCollection(release='latest'))
                ),
			    mode='minor', 
			    description='iljul99')

