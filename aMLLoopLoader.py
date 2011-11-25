"""

"""

import os, csv
from aMLSqlAlchemyClasses import *

# all loops
ifn = '/FR3D/MotifAtlas/Tables/All_loops.csv'

reader = csv.reader(open(ifn, 'rb'), delimiter=',', quotechar='"')
for row in reader:
    loop = AllLoops(parts = row)
    l = session.query(AllLoops).filter(AllLoops.id == loop.id).first()
    if not l:
        session.add(loop)

session.flush()
print 'Loops imported'

# loop modifications
ifn = '/FR3D/MotifAtlas/Tables/Loops_modifications.csv'
reader = csv.reader(open(ifn, 'rb'))
for row in reader:
    mod = LoopModifications(id=row[0], modification=row[1])
    m = session.query(LoopModifications).filter(LoopModifications.id == mod.id).first()
    if not m:
        session.add(mod)

session.flush()
print 'Loop modifications imported'

# loop qa
ifn = '/FR3D/MotifAtlas/Tables/Loops_QA.csv'
reader = csv.reader(open(ifn, 'rb'))

release = LoopRelease(mode='minor',description='test')
for row in reader:
    qa = LoopQA(id=row[0], code=row[1], release_id=release.id)
    q = session.query(LoopQA).filter(LoopQA.id == qa.id) \
                             .filter(LoopQA.release_id == release.id).first()
    if not q:
        session.add(qa)

#     session.merge(qa)

session.flush()
session.add(release)
session.commit()
print 'Loop QA imported'
