"""

A python script for managing the RNA 3D Hub cache.

Usage:
python manage_cache.py rna3dhub [option]
python manage_cache.py rna3dhub_dev [option]

where [option] is nrlist/motifs/loops/pdb

"""

import hashlib
import os
import urllib2
import sys
import logging

from sqlalchemy import distinct


from PdbInfoLoader import PdbInfoLoader
from MotifAtlasBaseClass import MotifAtlasBaseClass
from models import AllLoops, Motif, Release, NR_release, NR_class, session


class CacheManager(MotifAtlasBaseClass):

    def __init__(self, env):
        """
            env = rna3dhub or rna3dhub_dev
        """
        MotifAtlasBaseClass.__init__(self)
        self.cache_dir = os.path.join('/Servers', 'rna.bgsu.edu', env ,'application', 'cache')
        if not os.path.exists(self.cache_dir):
            os.mkdir(self.cache_dir)
        self.baseurl = 'http://rna.bgsu.edu/' + env

    def update_pdb_cache(self):
        """
           get a list of pdb files, loop over them, refresh cache
        """
        self._refresh_cache('pdb')
        P = PdbInfoLoader()
        P.get_all_rna_pdbs()
        for pdb in P.pdbs:
            logging.info(pdb)
            subpages = [pdb, pdb + '/motifs']
            for subpage in subpages:
                self._refresh_cache('pdb/' + subpage)

    def update_nrlist_cache(self):
        """
            get all nr_releases and nr_classes, then refresh cache
        """
        # homepage
        self._refresh_cache('nrlist')
        self._refresh_cache('nrlist/release/current')

        # http://rna.bgsu.edu/rna3dhub/nrlist/release/0.87/1.5A
        resolutions = ['1.5A', '2.0A', '2.5A', '3.0A', '3.5A', '4.0A', '20.0A', 'all']
        for release in session.query(NR_release).all():
            self._refresh_cache('nrlist/release/' + release.id)
            for resolution in resolutions:
                self._refresh_cache('nrlist/release/' + release.id + '/' + resolution)

        # http://rna.bgsu.edu/rna3dhub/nrlist/view/NR_4.0_81883.10
        for eq_class in session.query(distinct(NR_class.id)).all():
            self._refresh_cache('nrlist/view/' + eq_class[0])

    def update_loop_cache(self):
        """
            Delete old cache so that the pages can be regenerated with new data.
        """
        # http://rna.bgsu.edu/rna3dhub/loops/view/IL_157D_001
        for loop in session.query(AllLoops).all():
            self._refresh_cache('loops/view/' + loop.id, False)

    def update_motif_cache(self):
        """
            Refresh cache
        """
        # http://rna.bgsu.edu/rna3dhub/motifs/release/IL/current
        self._refresh_cache('motifs/release/IL/current')
        # http://rna.bgsu.edu/rna3dhub/motifs/release/HL/current
        self._refresh_cache('motifs/release/HL/current')

        # http://rna.bgsu.edu/rna3dhub/motifs/release/IL/0.6
        for release in session.query(Release).all():
            self._refresh_cache('motifs/release/' + release.type + '/' + release.id)

        # http://rna.bgsu.edu/rna3dhub/motif/view/IL_44742.1
        for motif in session.query(distinct(Motif.id)).all():
            self._refresh_cache('motif/view/' + motif[0])

    def _refresh_cache(self, uri_string, reload=True):
        """
            pdb/ + subpage
             NB! due to a bug in CodeIgniter, the urls are cached as
            'http://rna.bgsu.edu/rna3dhubpdb/1S72'
            no leading slash in uri_string
            example: pdb/1S72/motifs
            Reloads the page by default.
        """
        ci_url = self.baseurl + uri_string
        url = self.baseurl + '/' + uri_string
        md5 = hashlib.md5(ci_url).hexdigest()
        cached_file = os.path.join(self.cache_dir, md5)
        logging.info(url)
        if os.path.exists( cached_file ):
            logging.info("Deleting file %s" % cached_file)
            os.remove(cached_file)
        else:
            logging.info('File %s doesnt exist' % cached_file)
        if reload:
            try:
                response = urllib2.urlopen( url )
            except:
                logging.error('Failed to retrieve %s' % url)


def main(argv):
    """
    """

    C = CacheManager( argv[0] )
    C.start_logging()

    if argv[1] == 'pdb':
        C.update_pdb_cache()
    elif argv[1] == 'nrlist':
        C.update_nrlist_cache()
    elif argv[1] == 'motifs':
        C.update_motif_cache()
    elif argv[1] == 'loops':
        C.update_loop_cache()
    elif argv[1] == 'all':
        C.update_pdb_cache()
        C.update_nrlist_cache()
        C.update_motif_cache()
        C.update_loop_cache()
    else:
        logging.critical("Unrecognized option")
        sys.exit(1)

    C.set_email_subject('%s cache updated' % argv[1])
    C.send_report()


if __name__ == "__main__":
    main(sys.argv[1:])
