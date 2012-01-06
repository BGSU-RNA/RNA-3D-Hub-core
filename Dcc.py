import pdb, sys, getopt, os, re
from MLSqlAlchemyClasses import *
from sqlalchemy import distinct

# pdbfiles = ['1S72','2ZJR']

pdbfiles = set(session.query(func.substr(Dcc.id,1,4)).all())

# print pdbfiles

destination = '/Servers/rna.bgsu.edu/img/MotifAtlas/dcc_files/'


for pdbfile in pdbfiles:

    print pdbfile[0]
    nts = session.query(Coordinates).filter(Coordinates.nt_id.like(pdbfile[0]+'%')) \
                                    .order_by(Coordinates.nt_id).all()

    f = open(destination + pdbfile[0]+'.pdb', 'w')

    for nt in nts:
        parts = nt.nt_id.split('_')
        if len(parts[5]) == 1:
            f.write(nt.coordinates)
    f.close()

    f = open(destination + pdbfile[0]+'.dcc', 'w')

    l = dict()
    l['C'] = 24
    l['U'] = 23
    l['A'] = 26
    l['G'] = 27
    for coord, dcc in session.query(Coordinates, Dcc) \
                             .filter(Coordinates.nt_id == Dcc.id) \
                             .filter(Coordinates.nt_id.like(pdbfile[0]+'%')) \
                             .order_by(Coordinates.nt_id) \
                             .all():
        parts = dcc.id.split('_')

        for i in range(l[parts[5]]):
            f.write(str(dcc.sfcheck_correlation) + ' ' + \
                    str(dcc.sfcheck_correlation_side_chain) + ' ' + \
                    str(dcc.sfcheck_real_space_R) + ' ' + \
                    str(dcc.sfcheck_real_space_R_side_chain) + ' ' + \
                    str(dcc.sfcheck_connect) + ' ' + \
                    str(dcc.sfcheck_shift) + ' ' + \
                    str(dcc.sfcheck_shift_side_chain) + ' ' + \
                    str(dcc.sfcheck_density_index_main_chain) + ' ' + \
                    str(dcc.sfcheck_density_index_side_chain) + ' ' + \
                    str(dcc.sfcheck_B_iso_main_chain) + ' ' + \
                    str(dcc.sfcheck_B_iso_side_chain) + ' ' + \
                    str(dcc.mapman_correlation) + ' ' + \
                    str(dcc.mapman_real_space_R) + ' ' + \
                    str(dcc.mapman_Biso_mean) + ' ' + \
                    str(dcc.mapman_occupancy_mean) + '\n')


    f.close()

# zap; load auto "file:///FR3D/MotifAtlas/PDBDatabase/1S72/IL_1S72_028.pdb";if (_loadScript = '' && defaultLoadScript == '' && _filetype == 'Pdb') { select protein or nucleic;cartoons Only;color structure; select * }
# wireframe on
# cartoon off
# x='1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 3 3 3 3 3 3 3 3 3'
# DATA "property_dcc @x"
# color atoms property_dcc
# color property_dcc
# color {1.1} property_dcc
# color {1.2} property_dcc
# dots on
# dots off
# select [U]
# dots on
# center [U]
# center 1.1
# select {*}
# select {[U]}
# {[U]}.property_dcc=@x
# color [U] property_dcc
# color atoms property_dcc
# plot rama
# zap; load auto "file:///Servers/rna.bgsu.edu/nrlist/pdb/1S72.pdb";if (_loadScript = '' && defaultLoadScript == '' && _filetype == 'Pdb') { select protein or nucleic;cartoons Only;color structure; select * }
# select ({0:90984}) & connected; wireframe only;select ({90985:99038}) & !connected;stars 0.5;select *
# select protein;hide protein
# hide hetero
# hide protein
# select resno 2
# select 138
# navigate {138}
# navigate 2 center {138}
# cartoon off
# wireframe on
# hide hetero
# select 138
# select within(10, true, 138)
# hide *
# select within(10, true, 138)
# show
# display selected
# center 138
# zoom 100
# select rna
# display rna
# select 138.9
# select 138:9
# select 138:0
# navigate 2 {138:0}
# navigate 2 center {138:0}
# select within(10, true, 138:0)
# hide *
# display selected
# zoom 100
# center 138:0
# zoom 100
# zoom 300
# zoom 600
# zap
# zap; load auto "file:///Users/anton/Dropbox/Code/PyMotifLoader/temp.pdb";if (_loadScript = '' && defaultLoadScript == '' && _filetype == 'Pdb') { select protein or nucleic;cartoons Only;color structure; select * }
# select ({0:65154 65157:72508}) & connected; wireframe only;select ({65155 65156}) & !connected;stars 0.5;select *
# {*}.property_dcc=load('/Users/anton/Dropbox/Code/PyMotifLoader/dcc.dat')
# color atoms property_dcc
# color atoms property_dcc low
# color atoms property_dcc "low"
# color atoms property_dcc absolute 0.5 1.0
# color atoms property_dcc absolute 0.5 1.0 "low"
# color atoms property_dcc absolute 0.5 1.0
# color atoms property_dcc "low" absolute 0.5 1.0
# color atoms property_dcc "low" absolute 0.9 1.0
# color atoms property_dcc "low" absolute 0.95 1.0
# color atoms property_dcc "low" absolute 0.98 1.0
# zap
# zap; load auto "file:///Users/anton/Dropbox/Code/PyMotifLoader/2ZJR.pdb";if (_loadScript = '' && defaultLoadScript == '' && _filetype == 'Pdb') { select protein or nucleic;cartoons Only;color structure; select * }
# select ({0:70966}) & connected; wireframe only;select *
# {*}.property_dcc=data('/Users/anton/Code/PyMotifLoader/2ZJR.dcc')