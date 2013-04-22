from __future__ import with_statement
"""

Update pdb_info and pdb_obsolete tables. Pdb_info only contains the current
pdb files, old files are removed. Uses PDB RESTful services for getting a list
of all RNA-containing structures and for getting custom reports.

Usage: python PdbInfoLoader.py

"""

__author__ = 'Anton Petrov'

import csv
import logging
import traceback
import re
import urllib2
import sys
import os
from datetime import datetime
from ftplib import FTP
import pdb
from sqlalchemy import orm
import collections

from models import session, PdbInfo, PdbObsolete
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,parentdir)
from pdbx.reader.PdbxParser import PdbxReader
from MotifAtlasBaseClass import MotifAtlasBaseClass


class GetAllRnaPdbsError(Exception):
    """Raise when `get_all_rna_pdbs()` fails"""
    pass
class GetCustomReportError(Exception):
    """Raise when `_get_custom_report()` fails"""
    pass


class PdbInfoLoader(MotifAtlasBaseClass):
    """Class for retrieving RNA containing PDB files from the PDB."""

    def __init__(self):
        MotifAtlasBaseClass.__init__(self)
        self.pdbs = []
        self.adv_query_url     = 'http://www.rcsb.org/pdb/rest/search'
        self.custom_report_url = 'http://www.rcsb.org/pdb/rest/customReport'

    def update_pdb_info(self):
        """Get custom reports for all pdb files or for self.pdbs, if not empty,
           load into the database."""
        if not self.pdbs:
            self.get_all_rna_pdbs()
        failures = 0
        for pdb_id in self.pdbs:
            try:
                report = self._get_custom_report(pdb_id)
                self.__load_into_db(report)
            except:
                logging.error("Could not get and upate report for: %s", pdb_id)
                failures += 1
        if not failures:
            logging.info('Successful update of RNA-containing pdbs')
        else:
            logging.info('Partially successful update of RNA-containing pdbs')
            logging.info('Failed to update %s pdbs', failures)
        logging.info('%s', '+'*40)
        return True

    def get_all_rna_pdbs(self):
        """Get a list of all rna-containing pdb files, including hybrids. Raise
           a specific error if it fails."""
        logging.info('Getting a list of all rna-containing pdbs')
        query_text = """
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
        <containsProtein>I</containsProtein>
        <containsDna>I</containsDna>
        <containsRna>Y</containsRna>
        <containsHybrid>I</containsHybrid>
        </orgPdbQuery>
        """
        result = ''
        req = urllib2.Request(self.adv_query_url, data=query_text)
        f = urllib2.urlopen(req)
        result = f.read()

        if result:
            self.pdbs = result.split('\n')
            """filter out possible garbage in query results"""
            self.pdbs = filter(lambda x: len(x) == 4, self.pdbs)
            logging.info("Found %i PDB entries", len(self.pdbs))
        else:
            logging.critical("Failed to retrieve results")
            raise GetAllRnaPdbsError

    def _get_custom_report(self, pdb_id):
        """
            Gets a custom report in csv format for a single pdb file. Each chain
            is described in a separate line. Source organisms are retrived from
            cif files separately.
        """
        custom_report = '&customReportColumns=emResolution,structureId,chainId,structureTitle,experimentalTechnique,depositionDate,releaseDate,revisionDate,ndbId,resolution,classification,structureMolecularWeight,macromoleculeType,structureAuthor,entityId,sequence,chainLength,db_id,db_name,molecularWeight,secondaryStructure,entityMacromoleculeType,hetId,Ki,Kd,EC50,IC50,deltaG,deltaH,deltaS,Ka,compound,plasmid,taxonomyId,biologicalProcess,cellularComponent,molecularFunction,ecNo,expressionHost,cathId,cathDescription,scopId,scopDomain,scopFold,pfamAccession,pfamId,pfamDescription,crystallizationMethod,crystallizationTempK,phValue,densityMatthews,densityPercentSol,pdbxDetails,unitCellAngleAlpha,unitCellAngleBeta,unitCellAngleGamma,spaceGroup,lengthOfUnitCellLatticeA,lengthOfUnitCellLatticeB,lengthOfUnitCellLatticeC,Z_PDB,rObserved,rAll,rWork,rFree,refinementResolution,highResolutionLimit,reflectionsForRefinement,structureDeterminationMethod,conformerId,selectionCriteria,contents,solventSystem,ionicStrength,ph,pressure,pressureUnits,temperature,softwareAuthor,softwareName,version,method,details,conformerSelectionCriteria,totalConformersCalculated,totalConformersSubmitted,emdbId,emResolution,aggregationState,symmetryType,reconstructionMethod,specimenType&format=csv'
        url = '%s?pdbids=%s%s' % (self.custom_report_url, pdb_id, custom_report)
        logging.info('Getting custom report for %s', pdb_id)
        success = False
        retries = 0
        while success == False and retries < 3:
            try:
                f = urllib2.urlopen(url)
                result = f.read()
                success = True
            except:
                logging.warning(traceback.format_exc(sys.exc_info()))
                logging.warning("Failed to retrieve results")
                retries += 1
                result = None
                continue
        if result:
            logging.info("Retrieved custom report for %s", pdb_id)
        else:
            logging.critical("Failed to retrieve results")
            raise GetCustomReportError
        lines = result.split('<br />')
        return lines

    def _update_info(self, obj1, obj2):
        """obj1 - new instance as a dictionary, obj2 - old instance as a db object"""
        for k, v in obj1.iteritems():
            if k[0] == '_':
                continue # skip internal attributes
            setattr(obj2, k, v) # update the existing object
        return obj2

    def __load_into_db(self, lines):
        """Compares the custom report data with what's already in the database,
           reports any discrepancies and stores the most recent version."""
        logging.info('Loading custom report')
        description = lines.pop(0) # discard field names
        keys = description.split(',')

        for line in lines:
            """one line per chain"""
            if len(line) < 10:
                continue # skip empty lines
            """replace quotechars that are not preceded and followed by commas,
            except for the beginning and the end of the string
            example: in 1HXL unescaped doublequotes in the details field"""
            line = re.sub('(?<!^)(?<!,)"(?!,)(?!$)', "'", line)
            """parse the line using csv reader"""
            reader = csv.reader([line], delimiter=',', quotechar='"')
            """temporary store all field in a dictionary"""
            chain_dict = dict()
            for read in reader:
                for i, part in enumerate(read):
                    if not part:
                        part = None # to save as NULL in the db
                    chain_dict[keys[i]] = part
            """to save emResolution in resolution column"""
            if chain_dict['emResolution'] and not chain_dict['resolution']:
                chain_dict['resolution'] = chain_dict['emResolution']
            del(chain_dict['emResolution'])
            logging.info('%s %s', chain_dict['structureId'], chain_dict['chainId'])
            """check if this chain from this pdb is present in the db"""
            existing_chain = session.query(PdbInfo). \
                                     filter(PdbInfo.structureId==chain_dict['structureId']). \
                                     filter(PdbInfo.chainId==chain_dict['chainId']). \
                                     first()
            if existing_chain:
                """compare and merge objects"""
                existing_chain = self._update_info(chain_dict, existing_chain)
                session.merge(existing_chain)
            else:
                """create a new chain object"""
                new_chain = PdbInfo()
                for k, v in chain_dict.iteritems():
                    if not v:
                        v = None # to save as NULL in the db
                    setattr(new_chain, k, v)
                session.add(new_chain)

        session.commit()
        logging.info('Custom report saved in the database')

    def check_obsolete_structures(self):
        """Download the file with all obsolete structures over ftp, store
           the data in the database, remove obsolete entries
           from the pdb_info table"""

        TEMPFILE = 'obsolete.dat'
        REPEAT = 10
        done = False

        """download the data file from PDB"""
        for i in xrange(REPEAT):
            try:
                ftp = FTP('ftp.wwpdb.org')
                ftp.login()
                ftp.cwd('/pub/pdb/data/status')
                ftp.retrbinary("RETR %s" % TEMPFILE, open(TEMPFILE,"wb").write)
                ftp.quit()
                done = True
                logging.info('Downloaded obsolete.dat')
                break
            except Exception, e:
                logging.warning(e)
                logging.warning('Ftp download failed. Retrying...')

        if not done:
            logging.critical('All attempts to download obsolete.dat over ftp failed')
            logging.critical('Obsolete PDB files not updated')
            return

        """parse the data file"""
        obsolete_ids = []
        f = open(TEMPFILE, 'r')
        for line in f:
            # OBSLTE    26-SEP-06 2H33     2JM5 2OWI
            if 'OBSLTE' in line:
                parts = line.split()
                obsolete_date = datetime.strptime(parts[1], '%d-%b-%y')
                obsolete_ids.append(parts[2])
                replaced_by = ','.join(parts[3:])
                dbObj = PdbObsolete(obsolete_id=obsolete_ids[-1],
                                    date=obsolete_date,
                                    replaced_by=replaced_by)
                session.merge(dbObj)
        session.commit()

        """remove obsoleted files from pdb_info"""
        session.query(PdbInfo).\
                filter(PdbInfo.structureId.in_(obsolete_ids)).\
                delete(synchronize_session='fetch')

        """delete tempfile"""
        os.remove(TEMPFILE)

    def read_cif_file(self, pdb_id):
        """
            Open cif file with cif reader and return the object.
        """
        data = []
        filename = os.path.join(self.config['locations']['cif_dir'], pdb_id + '.cif')
        if not os.path.exists(filename):
            logging.info('Cif file %s not found' % filename)
            return None
        with open(filename, 'r') as raw:
            parser = PdbxReader(raw)
            parser.read(data)
        cif = data[0]
        return cif

    def _get_chain_id_map(self, cif):
        """
            Loop over the coordinates section to map chain ids to entity ids
            used internally in cif files. Inefficient, but guaranteed to work.
        """
        atom_site = cif.getObj('atom_site')
        chain_map = collections.defaultdict(dict)
        for i in xrange(atom_site.getRowCount()):
            chain     = atom_site.getValue('auth_asym_id', i)
            entity_id = atom_site.getValue('label_entity_id', i)
            if entity_id in chain_map and chain in chain_map[entity_id]:
                chain_map[entity_id][chain] += 1
            else:
                chain_map[entity_id][chain] = 1
        return chain_map

    def _get_source_organism_map(self, cif, pdb_id):
        """
            Map cif entity ids to organism names.
            entity_src_nat:
            Scientific name of the organism of the natural source.
            pdbx_entity_src_syn:
            Scientific name of the organism from which the entity was isolated.

            There are two additional fields related to source organisms:
            _pdbx_entity_src_syn.organism_scientific
            _em_entity_assembly.ebi_organism_scientific
            They seem to be empty for all RNA-containing 3D structures.
        """
        records = {'entity_src_nat':      'pdbx_organism_scientific', # Corresponds to SOURCE
                   'pdbx_entity_src_syn': 'organism_scientific',      # chemically synthesized
                   'entity_src_gen':      'pdbx_gene_src_scientific_name'}
        organism_map = dict()
        found = False
        for cif_category, cif_item in records.iteritems():
            block = cif.getObj(cif_category)
            if block is None:
                continue
            else:
                found = True
            for i in xrange(block.getRowCount()):
                organism  = block.getValue(cif_item, i)
                entity_id = block.getValue('entity_id', i)
                organism_map[entity_id] = organism
        if not found:
            logging.info('No cif source organisms for %s' % pdb_id)
        return organism_map

    def _get_organisms_by_chain(self, cif, pdb_id):
        """
            Map chain ids to organism names.
        """
        organisms = dict()
        organism_map = self._get_source_organism_map(cif, pdb_id)
        if not organism_map:
            return organisms # if no source organisms, return immediately
        chain_map = self._get_chain_id_map(cif)
        for entity_id, chains in chain_map.iteritems():
            for chain, val in chains.iteritems():
                if entity_id in organism_map:
                    organisms[chain] = organism_map[entity_id]
        pdb.set_trace()
        return organisms

    def __store_source_organisms(self, pdb_id, organisms):
        """
            Save source organisms in the database.
        """
        for chain, organism in organisms.iteritems():
            try:
                pdb_info = session.query(PdbInfo).\
                                   filter(PdbInfo.structureId==pdb_id).\
                                   filter(PdbInfo.chainId==chain).\
                                   one()
            except orm.exc.NoResultFound:
                """The pdb_info table is based on the information retrieved
                    over REST. For a small number of files not all chains are
                    reported, which causes this error. For example, 1ML5"""
                logging.info('PdbInfo does not have %s chain %s' %(pdb_id, chain))
                continue
            if organism == '?':
                organism = 'synthetic'
            pdb_info.source = organism
            session.merge(pdb_info)
        session.commit()

    def _get_organism_names_from_cif(self, pdb_id):
        """
            There are three ways of getting source organisms data:
            (1) parsing PDB files (2) using RESTful services (3) reading cif files.
            RESTful services don't return organisms for synthetic constructs.
            Parsing PDB files is unreliable, and using external tools
            introduces additional dependencies and may also be unreliable.
            Here we get the data from cif files using the pdbx parser that
            is also used for generating unit ids.
        """
        try:
            cif = self.read_cif_file(pdb_id)
            if cif:
                organisms = self._get_organisms_by_chain(cif, pdb_id)
                self.__store_source_organisms(pdb_id, organisms)
        except:
            logging.warning(traceback.format_exc(sys.exc_info()))
            logging.warning('There was a problem with %s' % pdb_id)

    def update_organism_names(self):
        """
            Convenience wrapper function.
        """
        logging.info('Updating organism names')
        for pdb_id in self.pdbs:
            logging.info(pdb_id)
            self._get_organism_names_from_cif(pdb_id)
        logging.info('Successful update of organism names')
        logging.info('%s', '+'*40)


def main(argv):
    """
    """
    logging.basicConfig(level=logging.DEBUG)

    P = PdbInfoLoader()

#     P.get_all_rna_pdbs()

    P.pdbs = ['1N8R']

#     P.update_pdb_info()
    P.update_organism_names()
#     P.check_obsolete_structures()


if __name__ == "__main__":
    main(sys.argv[1:])
