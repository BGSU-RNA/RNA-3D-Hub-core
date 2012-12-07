"""

This class stores information about non-redundant lists in a database.

"""

import logging
import sqlalchemy.exc


from models import session, NR_release, NR_class, NR_parents, NR_release_diff,\
                   NR_setdiff, NR_pdb, NR_handle
from MotifAtlasBaseClass import MotifAtlasBaseClass


class Uploader(MotifAtlasBaseClass):
    """
    When in release_diff mode, the class will not store the class ids, just the
    difference between releases
    release_mode: major/minor/reuse
    """
    def __init__(self, ensembles=None, release_mode='', release_description='', upload_mode=''):
        MotifAtlasBaseClass.__init__(self)
        self.c            = ensembles # collections, NRCollectionMerger
        self.motifs       = []
        self.loops        = []
        self.history      = []
        self.final_ids    = dict()
        self.intersection = []
        self.release_diff = []
        self.added_groups       = []
        self.removed_groups     = []
        self.updated_groups     = []
        self.old_updated_groups = []
        self.same_groups        = []
        self.added_pdbs         = []
        self.removed_pdbs       = []
        self.upload_mode = upload_mode
        self.release_mode = release_mode
        self.release_description = release_description

    def update_database(self):
        if self.upload_mode != 'release_diff':
            self.release = NR_release(mode=self.release_mode,
                                      description=self.release_description)
            self.__finalize()
            self.__process_set_diff()
            self.direct_parent = 1
            self.__process_release_diff()
            self.__commit()
        else:
            self.release = NR_release()
            self.release.id = self.c.c1.release
            self.direct_parent = 0
            self.__finalize()
            self.__process_release_diff()
            self.__release_diff_commit()

    def __finalize(self):
        """
        """
        self.added_pdbs   = self.c.c1.sl - self.c.c2.sl
        self.removed_pdbs = self.c.c2.sl - self.c.c1.sl

        for group_id in self.c.c1.sg:
            if group_id in self.c.new_ids:
                motif = NR_class(release_id=self.release.id,
                                 resolution=self.c.c1.res,
                                 comment=self.c.explanation[group_id])
                if self.c.parents.has_key(group_id):
                    parents = ','.join(set(self.c.parents[group_id]))
                else:
                    parents = ''
                if self.upload_mode == 'release_diff':
                    self.added_groups.append(group_id) # otherwise the new temporary id is committed
                else:
                    self.added_groups.append(motif.id)
                if self.upload_mode != 'release_diff':                    
	                logging.info('Group %s assigned new id %s' % (group_id, motif.id))

            elif group_id in self.c.correspond:
                old_id  = self.c.correspond[group_id]
                motif   = NR_class(id=old_id,
                                   release_id=self.release.id,
                                   increment=True,
                                   resolution=self.c.c1.res,
                                   comment=self.c.explanation[group_id])
                parents = ','.join(set([old_id] + self.c.parents[group_id]))
                self.updated_groups.append(motif.id)
                self.old_updated_groups.append(old_id)
                if self.upload_mode != 'release_diff':
                    logging.info('Group %s corresponds to motif %s and is assigned new id %s' % (group_id, old_id, motif.id))

            elif group_id in self.c.exact_match:
                id = self.c.exact_match[group_id]
                motif = NR_class(id=id,
                                 release_id=self.release.id,
                                 resolution=self.c.c1.res,
                                 comment=self.c.explanation[group_id])
                parents = ''
                self.same_groups.append(motif.id)
                # logging.info('Group %s matches exactly motif %s' % (group_id, motif.id))

            else:
                logging.error('Major problem')

            self.motifs.append(motif)
            self.final_ids[group_id] = motif.id

            if parents != '':
                self.history.append(NR_parents(class_id=motif.id, release_id=self.release.id, parents=parents))

            for loop_id in self.c.c1.d[group_id]:
                if loop_id in self.c.c1.reps:
                    self.loops.append(NR_pdb(id=loop_id, class_id=motif.id, release_id=self.release.id, rep = True))
                else:
                    self.loops.append(NR_pdb(id=loop_id, class_id=motif.id, release_id=self.release.id, rep = False))

            self.removed_groups = set(self.c.c2.groups) - set(self.old_updated_groups) - set(self.same_groups)

    def __process_release_diff(self):
        """
        """
        # do nothing during the very first upload
        if self.release.id == '0.1':
            return

        self.release_diff.append( NR_release_diff(
            nr_release_id1     = self.release.id,
            nr_release_id2     = self.c.c2.release,
            resolution         = self.c.c1.res,
            direct_parent      = self.direct_parent,
            added_groups       = ', '.join(self.added_groups),
            removed_groups     = ', '.join(self.removed_groups),
            updated_groups     = ', '.join(self.updated_groups),
            same_groups        = ', '.join(self.same_groups),
            added_pdbs         = ', '.join(self.added_pdbs),
            removed_pdbs       = ', '.join(self.removed_pdbs),
            num_added_groups   = len(self.added_groups),
            num_removed_groups = len(self.removed_groups),
            num_updated_groups = len(self.updated_groups),
            num_same_groups    = len(self.same_groups),
            num_added_pdbs     = len(self.added_pdbs),
            num_removed_pdbs   = len(self.removed_pdbs)
        ))


    def __process_set_diff(self):
        """
        """
        for loop_id in self.c.c1.sg:
            for motif_id in self.c.c2.sg:
                if motif_id != self.final_ids[loop_id] and \
                   self.c.intersection.has_key(loop_id) and \
                   self.c.intersection[loop_id].has_key(motif_id):

                    self.intersection.append( NR_setdiff(
                                nr_class1 = self.final_ids[loop_id],
                                nr_class2 = motif_id,
                                release_id = self.release.id,
                                intersection = ','.join(self.c.intersection[loop_id][motif_id]),
                                overlap = self.c.overlap[loop_id][motif_id],
                                one_minus_two = ','.join(self.c.setdiff[loop_id][motif_id]),
                                two_minus_one = ','.join(self.c.setdiff[motif_id][loop_id])
                    ))
                    self.intersection.append( NR_setdiff(
                                nr_class1 = motif_id,
                                nr_class2 = self.final_ids[loop_id],
                                release_id = self.release.id,
                                intersection = ','.join(self.c.intersection[loop_id][motif_id]),
                                overlap = self.c.overlap[motif_id][loop_id],
                                one_minus_two = ','.join(self.c.setdiff[motif_id][loop_id]),
                                two_minus_one = ','.join(self.c.setdiff[loop_id][motif_id])
                    ))

    def remove_release(self, release):
        """
        """
        try:
            session.query(NR_release).\
                    filter(NR_release.id==release).\
                    delete(synchronize_session='fetch')
            session.query(NR_class).\
                    filter(NR_class.release_id==release).\
                    delete(synchronize_session='fetch')
            session.query(NR_pdb).\
                    filter(NR_pdb.release_id==release).\
                    delete(synchronize_session='fetch')
            session.query(NR_setdiff).\
                    filter(NR_setdiff.release_id==release).\
                    delete(synchronize_session='fetch')
            session.query(NR_parents).\
                    filter(NR_parents.release_id==release).\
                    delete(synchronize_session='fetch')
            session.query(NR_release_diff).\
                    filter(NR_release_diff.nr_release_id1==release).\
                    delete(synchronize_session='fetch')
            session.query(NR_release_diff).\
                    filter(NR_release_diff.nr_release_id2==release).\
                    delete(synchronize_session='fetch')
            session.commit()
            logging.info('Release %s deleted successfully' % release)
        except:
            logging.error('Removing release %s failed' % release)
            session.rollback()
            sys.exit()

    def __release_diff_commit(self):
        """
        """
        try:
            session.add_all(self.release_diff)
            session.query(NR_handle).delete()
            session.commit()
            logging.info('Successful update')
        except sqlalchemy.exc.SQLAlchemyError, e:
            logging.error('Update failed. SQLAlchemy error. Rolling back.')
            logging.error(str(e))
            session.rollback()
            sys.exit()

    def __commit(self):
        """
        """
        try:
            handles = [i.handle for i in self.motifs]
            if len(handles) != len(set(handles)):
                pdb.set_trace()

            r = session.query(NR_release).filter(NR_release.id == self.release.id).first()
            if not r:
                session.add(self.release)

            session.add_all(self.motifs)
            session.add_all(self.loops)
            session.add_all(self.history)
            session.add_all(self.intersection)
            session.add_all(self.release_diff)

            session.commit()
            logging.info('Successful update')
        except sqlalchemy.exc.SQLAlchemyError, e:
            logging.error('Update failed. SQLAlchemy error. Rolling back.')
            logging.error(str(e))
            session.rollback()
            self.remove_release(self.release.id)
            sys.exit()
        except sqlalchemy.exc.DBAPIError, e:
            logging.error('Update failed. DBAPI error. Rolling back.')
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
