"""

Program for storing correspondences between the old and new style unit ids.

Example:
python UnitIdLoader.py

"""

__author__ = 'Anton Petrov'

import os
import csv
import pdb
import sys
import getopt
import logging


from models import session, PdbUnitIdCorrespondence, PdbAnalysisStatus
from MotifAtlasBaseClass import MotifAtlasBaseClass


# add dependencies to the python path
id_translation_submodule = 'UnitIdTranslation'
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.sys.path.append(parent_dir)
os.sys.path.append(os.path.join(parent_dir, id_translation_submodule))

# now import from the UnitIdTranslation folder
import id_translate as idt


class UnitIdLoader(MotifAtlasBaseClass):
    """
    """
    def __init__(self):
        MotifAtlasBaseClass.__init__(self)
        self.pdb_file_types   = ['.pdb', '.pdb1']
        self.pdb_files_folder = os.path.join(self.config['locations']['fr3d_root'],
                                             'FR3D',
                                             'PDBFiles')

    def import_unit_ids(self, pdbs, recalculate=False):
        """
            Determines what files need to be analyzed, deletes stored data if
            necessary, loops over the pdbs, generates unit ids and stores them
            in the database.
        """
        logging.info('Inside import_unit_ids')
        if not recalculate:
            recalculate = self.config['recalculate']['unit_ids']
        if recalculate:
            pdb_list = pdbs
            self.__delete_unit_ids(pdbs)
        else:
            pdb_list = self.filter_out_analyzed_pdbs(pdbs,'unit_ids')

        for pdb_id in pdb_list:
            for file_type in self.pdb_file_types:
                logging.info('Analyzing %s%s' % (pdb_id, file_type))
                pdb_file = os.path.join(self.pdb_files_folder, pdb_id + file_type)
                cif_file = os.path.join(self.pdb_files_folder, pdb_id + '.cif')

                skipped = False
                if not os.path.exists(pdb_file):
                    logging.warning('Skipping %s because %s%s was not found in %s' % (pdb_id, pdb_id, file_type, self.pdb_files_folder))
                    skipped = True
                    continue
                elif not os.path.exists(cif_file):
                    logging.warning('Skipping %s because %s.cif was not found in %s' % (pdb_id, pdb_id, self.pdb_files_folder))
                    skipped = True
                    continue

                try:
                    unit_ids = idt.get_id_correspondences(pdb_file, cif_file)
                except idt.LooksLikeAVirusStructureError:
                    logging.warning('%s looks like a viral structure' % pdb_id)
                    continue
                except Exception, err:
                    logging.critical('Crash on %s' % pdb_id)
                    continue
                logging.info('Found %i id pairs' % len(unit_ids))

                self.__store_unit_ids(unit_ids)
            if not skipped:
                self.mark_pdb_as_analyzed(pdb_id, 'unit_ids')

        logging.info('%s', '='*40)

    def __delete_unit_ids(self, pdbs):
        """
        """
        logging.info('Deleting unit id correspondence data')
        session.query(PdbUnitIdCorrespondence). \
                filter(PdbUnitIdCorrespondence.pdb.in_(pdbs)). \
                delete(synchronize_session='fetch')
        logging.info('Resetting unit id analysis status')
        for statusObj in session.query(PdbAnalysisStatus). \
                                 filter(PdbAnalysisStatus.id.in_(pdbs)). \
                                 all():
            statusObj.unit_ids = None
            session.merge(statusObj)
        session.commit()

    def __store_unit_ids(self, unit_ids):
        """
        """
        logging.info('Importing unit ids')
        for (old_id, new_id) in unit_ids:
            # old id 1EKA_AU_1_B_8_C_
            (pdb_id,au_ba,model,chain,seq_id,comp_id,ins_code) = old_id.split('_')
            pdb_file = 'pdb' if au_ba == 'AU' else 'pdb' + au_ba[2:] # BA1, BA10

            underscores = new_id.count('_')
            if underscores == 8: # all fields present
                (pdb_id,model,chain,comp_id,seq_id,alt_id,ins_code,sym_op1,sym_op2) = new_id.split('_')
                sym_op = sym_op1 + '_' + sym_op2
            elif underscores == 6: # default sym_op
                (pdb_id,model,chain,comp_id,seq_id,alt_id,ins_code) = new_id.split('_')
                sym_op = '1_555'
            elif underscores == 5: # default sym_op and ins_code
                (pdb_id,model,chain,comp_id,seq_id,alt_id) = new_id.split('_')
                sym_op   = '1_555'
                ins_code = ''
            elif underscores == 4: # default sym_op, ins_code, alt_id
                (pdb_id,model,chain,comp_id,seq_id) = new_id.split('_')
                sym_op   = '1_555'
                ins_code = ''
                alt_id   = ''
            else:
                msg = 'Unknown id format %s' % new_id
                logging.critical(msg)
                session.rollback()
                self._crash(msg)

            U = PdbUnitIdCorrespondence(old_id   = old_id,
                                        unit_id  = new_id,
                                        pdb      = pdb_id,
                                        model    = model,
                                        chain    = chain,
                                        seq_id   = seq_id,
                                        comp_id  = comp_id,
                                        alt_id   = alt_id,
                                        ins_code = ins_code,
                                        sym_op   = sym_op,
                                        pdb_file = pdb_file)
            try:
                session.add(U)
            except:
                pass
        session.commit()
        logging.info('Ids successfully imported')



def usage():
    print __doc__

def main(argv):
    """
    """
    U = UnitIdLoader()
    U.start_logging()

    from PdbInfoLoader import PdbInfoLoader
    P = PdbInfoLoader()
    P.get_all_rna_pdbs()

    try:
        U.import_unit_ids(P.pdbs)
    except:
        e = sys.exc_info()[1]
        U.set_email_subject('Unit id update failed')
        U._crash(e)

    U.set_email_subject('Unit ids successfully updated')
    U.send_report()

if __name__ == "__main__":
    main(sys.argv[1:])