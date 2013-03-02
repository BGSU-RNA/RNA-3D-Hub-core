"""


"""

import logging
import csv
import os
import shutil
import re
import sys

from sqlalchemy.sql import or_
from sqlalchemy import *
import sqlalchemy.exc

from MotifAtlasBaseClass import MotifAtlasBaseClass
from models import Release, session, ML_handle, LoopOrder, LoopPosition, \
                   Motif, Loop, SetDiff, Parents, LoopDiscrepancy, Release_diff, \
                   MotifAnnotation



class Uploader(MotifAtlasBaseClass):
    """
    """
    def __init__(self, ensembles=None, mode='', description='', files={}, upload_mode='', motif_type=''):
        self.c           = ensembles
        self.motif_type  = motif_type.upper()
        self.release     = Release(mode=mode, description=description, type=self.motif_type)
        self.upload_mode = upload_mode
        self.motifs    = []
        self.loops     = []
        self.history   = []
        self.intersection = []
        self.final_ids    = dict()
        self.files        = files

        self.release_diff       = []
        self.added_groups       = []
        self.removed_groups     = []
        self.updated_groups     = []
        self.old_updated_groups = []
        self.same_groups        = []
        self.added_loops        = []
        self.removed_loops      = []
        self.loop_order       = []
        self.loop_positions   = []
        self.loop_discrepancy = []
        self.bp_signatures = []

    def import_release(self):
        """
        """
        if self.upload_mode != 'release_diff':
            self.direct_parent = 1
            self.__finalize()
            self.__process_set_diff()
            self.__process_release_diff()
            self.__process_motif_loop_order()
            self.__process_motif_loop_positions()
            self.__process_mutual_discrepancy()
            self.__process_bp_signatures()
            self.__store_correspondences()
            self.__commit()
            self.__rename_mat_files()
            self.__move_2d_files()
            self.__process_graphml_file()
        else:
            self.direct_parent = 0
            self.release = Release()
            self.release.id = self.c.c1.release
            self.__finalize()
            self.__process_release_diff()
            self.__release_diff_commit()

    def __finalize(self):
        """
        """
        self.added_loops   = self.c.c1.sl - self.c.c2.sl
        self.removed_loops = self.c.c2.sl - self.c.c1.sl

        for group_id in self.c.c1.sg:
            if group_id in self.c.new_ids:
                motif = Motif(release_id=self.release.id,
                              comment=self.c.explanation[group_id],
                              type=self.motif_type)
                if self.c.parents.has_key(group_id):
                    parents = ','.join(set(self.c.parents[group_id]))
                else:
                    parents = ''
                if self.upload_mode == 'release_diff':
                    self.added_groups.append(group_id) # otherwise the new temporary id is committed
                else:
                    self.added_groups.append(motif.id)
                logging.info('Group %s assigned new id %s' % (group_id, motif.id))
            elif group_id in self.c.correspond:
                old_id  = self.c.correspond[group_id]
                motif   = Motif(id=old_id,
                                release_id=self.release.id,
                                increment=True,
                                comment=self.c.explanation[group_id],
                                type=self.motif_type)
                parents = ','.join(set([old_id] + self.c.parents[group_id]))
                self.updated_groups.append(motif.id)
                self.old_updated_groups.append(old_id)
                logging.info('Group %s corresponds to motif %s and is assigned new id %s' % (group_id, old_id, motif.id))
            elif group_id in self.c.exact_match:
                id = self.c.exact_match[group_id]
                motif = Motif(id=id,
                              release_id=self.release.id,
                              comment=self.c.explanation[group_id],
                              type=self.motif_type)
                self.same_groups.append(motif.id)
                parents = ''
                logging.info('Group %s matches exactly motif %s' % (group_id, motif.id))
            else:
                logging.error('Major problem')

            self.motifs.append(motif)
            if parents != '':
                if self.upload_mode != 'release_diff':
                    self.__inherit_history(motif, parents)
                self.history.append(Parents(motif_id=motif.id, release_id=self.release.id, parents=parents))
            self.final_ids[group_id] = motif.id
            for loop_id in self.c.c1.d[group_id]:
                self.loops.append(Loop(id=loop_id, motif_id=motif.id, release_id=self.release.id))

            self.removed_groups = set(self.c.c2.groups) - set(self.old_updated_groups) - set(self.same_groups)

    def __store_correspondences(self):
        """
        """
        f = open(self.files['correspondences'], 'w')
        for group_id, motif_id in self.final_ids.iteritems():
            f.write('%s,%s\n' % (group_id, motif_id) )
        f.close()

    def __inherit_history(self, motif, parents):
        """
        """
        parent_motifs = parents.split(',')
        common_name = []
        annotation = []
        author = []
        for parent in parent_motifs:
            parent_motif = session.query(MotifAnnotation).\
                                   filter(MotifAnnotation.motif_id==parent).\
                                   first()
            if parent_motif:
                if parent_motif.common_name is not None and parent_motif.common_name != '':
                    common_name.append(parent_motif.common_name)
                if parent_motif.annotation is not None and parent_motif.annotation != '':
                    annotation.append(parent_motif.annotation)
                if parent_motif.author is not None and parent_motif.author != '':
                    author.append(parent_motif.author)

        session.merge(MotifAnnotation(motif_id    = motif.id,
                                      common_name = ' | '.join(set(common_name)),
                                      annotation  = ' | '.join(set(annotation)),
                                      author      = ' | '.join(set(author))))

    def __process_release_diff(self):
        """
        """
        # do nothing during the very first upload
        if self.release.id == '0.1':
            return

        self.release_diff.append( Release_diff(
            release_id1        = self.release.id,
            release_id2        = self.c.c2.release,
            type               = self.motif_type,
            direct_parent      = self.direct_parent,
            added_groups       = ', '.join(sorted(self.added_groups)),
            removed_groups     = ', '.join(sorted(self.removed_groups)),
            updated_groups     = ', '.join(sorted(self.updated_groups)),
            same_groups        = ', '.join(sorted(self.same_groups)),
            added_loops        = ', '.join(sorted(self.added_loops)),
            removed_loops      = ', '.join(sorted(self.removed_loops)),
            num_added_groups   = len(self.added_groups),
            num_removed_groups = len(self.removed_groups),
            num_updated_groups = len(self.updated_groups),
            num_same_groups    = len(self.same_groups),
            num_added_loops    = len(self.added_loops),
            num_removed_loops  = len(self.removed_loops)
        ))

    def __process_set_diff(self):
        """
        """
        for loop_id in self.c.c1.sg:
            for motif_id in self.c.c2.sg:
                if motif_id != self.final_ids[loop_id] and \
                   self.c.intersection.has_key(loop_id) and \
                   self.c.intersection[loop_id].has_key(motif_id):

                    self.intersection.append( SetDiff(
                                motif_id1 = self.final_ids[loop_id],
                                motif_id2 = motif_id,
                                release_id = self.release.id,
                                intersection = ','.join(self.c.intersection[loop_id][motif_id]),
                                overlap = self.c.overlap[loop_id][motif_id],
                                one_minus_two = ','.join(self.c.setdiff[loop_id][motif_id]),
                                two_minus_one = ','.join(self.c.setdiff[motif_id][loop_id])
                    ))
                    self.intersection.append( SetDiff(
                                motif_id1 = motif_id,
                                motif_id2 = self.final_ids[loop_id],
                                release_id = self.release.id,
                                intersection = ','.join(self.c.intersection[loop_id][motif_id]),
                                overlap = self.c.overlap[motif_id][loop_id],
                                one_minus_two = ','.join(self.c.setdiff[motif_id][loop_id]),
                                two_minus_one = ','.join(self.c.setdiff[loop_id][motif_id])
                    ))

    def remove_release(self, release):
        """
           Must pay attention to release type so that only IL or HL release
           is deleted, not both
        """
        # ml_handles
        logging.info('Removing release %s from table %s' % (release, ML_handle.__tablename__) )
        session.query(ML_handle).delete(synchronize_session='fetch')
        # ml_loop_order
        logging.info('Removing release %s from table %s' % (release, LoopOrder.__tablename__) )
        session.query(LoopOrder).\
                filter(LoopOrder.release_id==release).\
                filter(LoopOrder.motif_id.like(self.motif_type+'%')).\
                delete(synchronize_session='fetch')
        # ml_loop_positions
        logging.info('Removing release %s from table %s' % (release, LoopPosition.__tablename__) )
        session.query(LoopPosition).\
                filter(LoopPosition.release_id==release).\
                filter(LoopPosition.motif_id.like(self.motif_type+'%')).\
                delete(synchronize_session='fetch')
        # ml_releases
        logging.info('Removing release %s from table %s' % (release, Release.__tablename__) )
        session.query(Release).\
                filter_by(id=release).\
                filter_by(type=self.motif_type).\
                delete(synchronize_session='fetch')
        # ml_motifs
        logging.info('Removing release %s from table %s' % (release, Motif.__tablename__) )
        session.query(Motif).\
                filter_by(release_id=release).\
                filter_by(type=self.motif_type).\
                delete(synchronize_session='fetch')
        # ml_loops
        logging.info('Removing release %s from table %s' % (release, Loop.__tablename__) )
        session.query(Loop).\
                filter(Loop.release_id==release).\
                filter(Loop.id.like(self.motif_type+'%')).\
                delete(synchronize_session='fetch')
        # ml_set_diff
        logging.info('Removing release %s from table %s' % (release, SetDiff.__tablename__) )
        session.query(SetDiff).\
                filter(SetDiff.release_id==release).\
                filter(SetDiff.motif_id1.like(self.motif_type+'%')).\
                delete(synchronize_session='fetch')
        # ml_history
        logging.info('Removing release %s from table %s' % (release, Parents.__tablename__) )
        session.query(Parents).\
                filter(Parents.release_id==release).\
                filter(Parents.motif_id.like(self.motif_type+'%')).\
                delete(synchronize_session='fetch')
        # ml_mutual_discrepancy
        logging.info('Removing release %s from table %s' % (release, LoopDiscrepancy.__tablename__) )
        session.query(LoopDiscrepancy).\
                filter(LoopDiscrepancy.release_id==release).\
                delete(synchronize_session='fetch')
        # ml_release_diff
        logging.info('Removing release %s from table %s' % (release, Release_diff.__tablename__) )
        session.query(Release_diff).\
                filter(or_(Release_diff.release_id1==release, Release_diff.release_id2==release)).\
                filter_by(type=self.motif_type).\
                delete(synchronize_session='fetch')
        session.commit()

    def __release_diff_commit(self):
        """
        """
        try:
            session.add_all(self.release_diff)
            session.query(ML_handle).delete()
            session.commit()
            logging.info('Successful update')
        except sqlalchemy.exc.SQLAlchemyError, e:
            logging.error('Update failed. SQLAlchemy error. Rolling back.')
            logging.error(str(e))
            session.rollback()
            sys.exit()

    def __commit(self):
        try:
            session.add(self.release)
            session.add_all(self.motifs)
            session.add_all(self.loops)
            session.add_all(self.history)
            session.add_all(self.intersection)
            session.add_all(self.release_diff)
            session.add_all(self.loop_order)
            session.add_all(self.loop_positions)
            session.add_all(self.loop_discrepancy)
            session.commit()
            logging.info('Successful update')
        except sqlalchemy.exc.SQLAlchemyError, e:
            logging.error('Update failed. Rolling back.')
            logging.error(str(e))
            session.rollback()
            self.remove_release(self.release.id)
            sys.exit()
        except sqlalchemy.exc.DBAPIError, e:
            logging.error('Update failed. Rolling back.')
            logging.error(str(e))
            session.rollback()
            self.remove_release(self.release.id)
            sys.exit()
        except sys.exc_info()[0]:
            logging.error('Update failed. Rolling back.')
            logging.error(sys.exc_info()[0])
            session.rollback()
            self.remove_release(self.release.id)
            sys.exit()

    def __process_motif_loop_order(self):
        """
        IL_017,IL_1O9M_003,1,1
        """
        r = csv.reader(open(self.files['MotifLoopOrder']), delimiter=',', quotechar='"')
        for row in r:
            self.loop_order.append(LoopOrder(motif_id=self.final_ids[row[0]],
                                        loop_id=row[1],
                                        release_id=self.release.id,
                                        original_order=row[2],
                                        similarity_order=row[3]
                                   ))


    def __process_motif_loop_positions(self):
        """
            IL_017,IL_1O9M_003,1O9M_AU_1_A_1491_G_,1
        """
        r = csv.reader(open(self.files['MotifPositions']), delimiter=',', quotechar='"')
        for row in r:
            self.loop_positions.append(LoopPosition(motif_id=self.final_ids[row[0]],
                                        loop_id=row[1],
                                        release_id=self.release.id,
                                        nt_id=row[2],
                                        position=row[3]
                                   ))


    def __process_bp_signatures(self):
        """
            "Group_001","06_cWW-cWW-cWW"
        """
        r = csv.reader(open(self.files['BpSignatures']), delimiter=',', quotechar='"')
        for row in r:
            session.merge(MotifAnnotation(motif_id=self.final_ids[row[0]],
                                          bp_signature=row[1]))
        session.commit()

    def __process_mutual_discrepancy(self):
        """
            IL_1O9M_003,0.0000,IL_1O9M_003
        """
        r = csv.reader(open(self.files['MutualDiscrepancy']), delimiter=',', quotechar='"')
        for row in r:
            self.loop_discrepancy.append(LoopDiscrepancy(loop_id1=row[0],
                                        discrepancy=row[1],
                                        loop_id2=row[2],
                                        release_id=self.release.id,
                                   ))

    def __rename_mat_files(self):
        """
        """
        if not os.path.exists(self.files['MatFiles_dest']):
            os.mkdir(self.files['MatFiles_dest'])

        for file in self.c.c1.sg:
            src = os.path.join(self.files['MatFiles_origin'], file+'.mat')
            dst = os.path.join(self.files['MatFiles_dest'], self.final_ids[file] + '.mat')
            if os.path.exists(src):
                shutil.copyfile(src, dst)
            else:
                logging.warning("File %s wasn't found" % src)

    def __move_2d_files(self):
        """
        """
        img_release_folder = os.path.join(self.files['2ds_destination'],self.motif_type + self.release.id)
        if not os.path.exists(img_release_folder):
            os.mkdir(img_release_folder)

        for file in self.c.c1.sg:
            src = os.path.join(self.files['2ds_origin'], file + '.png')
            dst = os.path.join(img_release_folder, self.final_ids[file] + '.png')
            if os.path.exists(src):
                shutil.copyfile(src, dst)
            else:
                logging.warning("File %s wasn't found" % src)

    def __process_graphml_file(self):
        """
        """
        graphml = os.path.join(self.files['folder'],'Supergroups.graphml')
        if os.path.exists(graphml):
            f = open(graphml, 'r')
            contents = f.read()
            for group in self.c.c1.sg:
#                 if group in contents:
#                     # first try replacing group names as they were given
#                     contents = contents.replace(group, self.final_ids[group])
#                 else:
#                     # try to replace ids like Group_001 in case the matfiles were renamed to IL_001 etc
                parts = group.split('_')
                contents = contents.replace('Group_'+parts[1], self.final_ids[group])
            contents = contents.replace('\n','')
            self.release.graphml = contents
            session.add(self.release)
            session.commit()
