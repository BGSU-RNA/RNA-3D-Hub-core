
from models import LoopSearch, session

f = open('/Users/anton/FR3D/loops.txt')

lines = f.readlines()

loops = lines[0].split(',')

N = len(loops)
todo = 0

for i, loop in enumerate(loops):
    print i
    todo += N - len(session.query(LoopSearch).filter_by(loop_id1=loop).filter(LoopSearch.loop_id2.in_(loops)).all())
    print todo

print todo
