"""
Read .cif files to get assembly information.
Read the assembly_gen and struct_oper tables
"""

# import requests
from pymotifs import core
from pymotifs import models as mod
# from pymotifs.utils.pdb import CustomReportHelper
from pymotifs.download import Downloader
from pymotifs.pdbs.loader import Loader as PdbLoader

from fr3d.cif.reader import Cif


class Loader(core.SimpleLoader):
    """
    """

    merge_data = True
    allow_no_data = True

    dependencies = set([Downloader, PdbLoader])

    @property
    def table(self):
        return mod.AssemblyInfo


    def to_process(self, pdbs, **kwargs):
        """
        Return all PDB ids the first time.

        :param list pdbs: The list of pdb ids of interest on this run.
        :param dict kwargs: The keyword arguments which are ignored.
        :returns: A list of pdb ids to process.
        """

        if len(pdbs) == 1:
            return pdbs

        # process all pdb files
        with self.session() as session:
            query = session.query(mod.PdbInfo.pdb_id).\
            distinct()
        pdbs = [row.pdb_id for row in query]

        return pdbs


    def has_data(self, pdb, **kwargs):
        """
        Check if we have data for the given PDB.

        Parameters
        ----------
        pdb : str
            The PDB to check

        Returns
        -------
        has_data : bool
            True if we have data for the PDB.
        """

        with self.session() as session:
            query = session.query(mod.AssemblyInfo).\
                filter_by(pdb_id=pdb).\
                limit(1)
            return bool(query.count())


    def query(self, session, pdb):
        """
        Generate a query to find all entries in PDBInfo for the given
        PDB id.  Added in November 2020 so this can be a SimpleLoader.
        Attributes
        ----------
        session : Session
            The `Session` to use.
        pdb : str
            The PDB id to use.
        Returns
        -------
        query : Query
            Returns an SqlAlchemy query for all entires in pdb_info with
            the given PDB id.
        """
        return session.query(mod.AssemblyInfo).\
            filter_by(pdb_id=pdb)


    def data(self, pdb, **kwargs):
        """
        Read the .cif file, extract the data we need, return database entries
        """

        # open entire cif file, not just residues
        cif_access = self.cif(pdb)

        # To get the values of oper_expression, read this table:
        symmetry_id_to_name = {}
        operators = cif_access.table('pdbx_struct_oper_list').rows
        for operator in operators:
            print(operator)
            id = operator['id']
            name = operator['name']
            symmetry_id_to_name[id] = name

        # t.rows is a list of dictionaries with these keys: 'assembly_id', 'oper_expression', 'asym_id_list'
        assemblies = cif_access.table('pdbx_struct_assembly_gen').rows

        for assembly in assemblies:
            print(assembly)

            chains = assembly["asym_id_list"].split(",")
            symmetry_ids = assembly["oper_expression"].split(",")

            for chain in chains:
                for symmetry_id in symmetry_ids:
                    data = {}
                    data["pdb_id"] = pdb
                    data["assembly_id"] = int(assembly['assembly_id'])     # integer
                    data["chain_name"] = chain
                    data["symmetry_id"] = int(symmetry_id)     # integer
                    data["symmetry"] = symmetry_id_to_name[symmetry_id]

                    print(data)

                    yield(mod.AssemblyInfo(**data))
