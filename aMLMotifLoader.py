"""

"""

import os
import aMLCollection as motifs
from aMLSqlAlchemyClasses import *


f = dict()
f['folder']            = '/FR3D/Workspace/Releases/iljul2'
f['MotifList']         = os.path.join(f['folder'],'MotifList.csv')
f['MotifLoopOrder']    = os.path.join(f['folder'],'MotifLoopOrder.csv')
f['MotifPositions']    = os.path.join(f['folder'],'MotifPositions.csv')
f['MutualDiscrepancy'] = os.path.join(f['folder'],'MutualDiscrepancy.csv')
f['MatFiles']          = os.path.join(f['folder'],'mat')

    
A = Uploader(collections=motifs.MotifCollectionMerger(
                    motifs.MotifCollection(file=f['MotifList']),
#                     motifs.MotifCollection(file='input/iljul4.csv'),                    
                    motifs.MotifCollection(collection=LoopCollection(release='latest'))
                ),
			    mode='minor', 
			    description='iljul4.csv',
			    files=f)

