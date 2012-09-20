"""

A python script for managing the RNA 3D Hub cache.

Usage:
python manage_cache.py rna3dhub
python manage_cache.py rna3dhub_dev

"""

import hashlib
import os
import urllib2
import sys
import logging

from PdbInfoLoader import PdbInfoLoader
from NRSqlAlchemyClasses import session, NR_release, NR_class
from MLSqlAlchemyClasses import AllLoops, Motif, Release

from sqlalchemy import distinct


class CacheManager():

    def __init__(self, env):
        """
            env = rna3dhub or rna3dhub_dev
        """
        self.cache_dir = os.path.join('/Servers', 'rna.bgsu.edu', env ,'application', 'cache')
        self.baseurl = 'http://rna.bgsu.edu/' + env

    def update_pdb_cache(self, pdbs):
        """
           get a list of pdb files, loop over them, delete existing cache and
           request a new version
        """
        for pdb in pdbs:
            logging.info(pdb)
            subpages = [pdb, pdb + '/motifs']
            for subpage in subpages:
                self._refresh_cache('pdb/' + subpage)

    def update_nrlist_cache(self):
        """
            get all nr_releases and nr_classes. Loop over and refresh cache
        """
        # homepage
        self._refresh_cache('nrlist')

        # http://rna.bgsu.edu/rna3dhub/nrlist/release/0.87
        for release in session.query(NR_release).all():
            self._refresh_cache('nrlist/release/' + release.id)

        # http://rna.bgsu.edu/rna3dhub/nrlist/view/NR_4.0_81883.10
        for eq_class in session.query(distinct(NR_class.id)).all():
            self._refresh_cache('nrlist/view/' + eq_class[0])

    def update_loop_cache(self):
        """
        """
        # http://rna.bgsu.edu/rna3dhub/loops/view/IL_157D_001
        for loop in session.query(AllLoops).all():
            self._refresh_cache('loops/view/' + loop.id)

    def update_motifs_cache(self):
        """
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
            self._refresh_cache('motifs/view' + motif.id)

    def _refresh_cache(self, uri_string):
        """
            pdb/ + subpage
             NB! due to a bug in CodeIgniter, the urls are cached as
            'http://rna.bgsu.edu/rna3dhubpdb/1S72'
            no leading slash in uri_string
            example: pdb/1S72/motifs
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
        try:
            response = urllib2.urlopen( url )
        except:
            logging.error('Failed to retrieve %s' % url)


def main(argv):
    """
    """
    logging.basicConfig(level=logging.DEBUG)

    P = PdbInfoLoader()
    P.get_all_rna_pdbs()

    C = CacheManager( argv[0] )

    C.update_pdb_cache( P.pdbs )
    C.update_nrlist_cache()
    C.update_motif_cache()
    C.update_loop_cache()


if __name__ == "__main__":
    main(sys.argv[1:])
