"""

"""

import os, csv
from aMLSqlAlchemyClasses import *

# all loops
ifn = '/FR3D/MotifAtlas/Tables/All_loops.csv'

reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
for row in reader:
    loop = AllLoops(parts = row)
    session.merge(loop)

session.flush()
print 'Loops imported'

# loop modifications
ifn = '/FR3D/MotifAtlas/Tables/Loops_modifications.csv'
reader = csv.reader(open(ifn, 'rb'))
for row in reader:
    mod = LoopModifications(id=row[0], modification=row[1])
    session.merge(mod)

session.flush()
print 'Loop modifications imported'

# loop qa
ifn = '/FR3D/MotifAtlas/Tables/Loops_QA.csv'
reader = csv.reader(open(ifn, 'rb'))

release = LoopRelease(mode='minor',description='test')
for row in reader:
    qa = LoopQA(id=row[0], code=row[1], release_id=release.id)
    session.merge(qa)

session.flush()
session.add(release)
session.commit()
print 'Loop QA imported'
