#!/usr/bin/python

import shutil, os, re
from MLSqlAlchemyClasses import session, MotifAnnotation


f = open('/Users/anton/FR3D/MotifAtlas/Releases/ilmarch/correspondences.txt', 'r')

from mlabwrap import mlab

for line in f:
        m = re.search(r"(Group_\d+).+(IL_\d+\.\d+)", line)
        if m:
            print "%s %s" % (m.group(1), m.group(2))
            s = mlab.getBasepairSignature('/Users/anton/FR3D/MotifAtlas/Releases/ilmarch/mat/' + m.group(2) + '.mat')
            print s
            m = MotifAnnotation(motif_id = m.group(2),
                                bp_signature = s)
            session.merge(m)
        session.commit()
#             src = os.path.join(os.getcwd(), 'Groups', m.group(1) + '.mat')
#             dst = os.path.join(os.getcwd(), 'mat', m.group(2) + '.mat')
#             shutil.copyfile(src, dst)