"""
Read .cif files to get assembly information.
Read the assembly_gen and struct_oper tables
"""

from pymotifs import core
from pymotifs import models as mod
from pymotifs.download import Downloader
from pymotifs.pdbs.info import Loader as InfoLoader

from fr3d.cif.reader import Cif


class Loader(core.SimpleLoader):

    merge_data = True
    allow_no_data = True

    dependencies = set([Downloader, InfoLoader])

    @property
    def table(self):
        return mod.AssemblyInfo


    def to_process(self, pdbs, **kwargs):
        """
        From among pdbs, find the ones that do not have assembly information
        in the assembly_info table.

        :param list pdbs: The list of pdb ids of interest on this run.
        :param dict kwargs: The keyword arguments, which are ignored.
        :returns: A list of pdb ids to process.
        """

        if len(pdbs) == 1:
            return pdbs

        # query to find pdb ids that already have assembly information
        with self.session() as session:
            query = session.query(mod.AssemblyInfo.pdb_id).\
                distinct()
            pdbs_with_assembly = set([row.pdb_id for row in query])

        pdbs_needed = set(pdbs) - pdbs_with_assembly

        if len(pdbs_needed) == 0:
            raise core.Skip("No pdbs need assembly data loaded")

        return sorted(pdbs_needed)


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
        PDB id.
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

        # loop over residues to build the mapping from PDB chain to author chain
        pdb_chain_to_author_chain = {}
        c = 0
        try:
            for atom_site in cif_access.table('atom_site'):
                pdb_chain = atom_site['label_asym_id']
                author_chain = atom_site['auth_asym_id']
                pdb_chain_to_author_chain[pdb_chain] = author_chain
                c += 1
        except:
            # some structures, like 9A8W, do not have an atom_site block,
            # but instead have ihm_starting_model_coord
            pass

        # print(sorted(pdb_chain_to_author_chain.items()))
        self.logger.info(sorted(pdb_chain_to_author_chain.items()))
        # print('Processed %d lines of ATOM records' % c)

        # To get the values of oper_expression, read this table
        # Unfortunately, some structures do not have this block.
        symmetry_id_to_name = {}
        assemblies = []
        try:
            operators = cif_access.table('pdbx_struct_oper_list').rows
            for operator in operators:
                self.logger.info(operator)
                id = operator['id']
                name = operator.get('name',None)
                symmetry_id_to_name[id] = name

            # t.rows is a list of dictionaries with these keys: 'assembly_id', 'oper_expression', 'asym_id_list'
            assemblies = cif_access.table('pdbx_struct_assembly_gen').rows

            max_assembly_id = 0
            for assembly in assemblies:
                # print("Working on assembly %s in %s" % (assembly,pdb))

                pdb_chains = assembly["asym_id_list"].split(",")

                try:
                    assembly_id = int(assembly['assembly_id'])
                except:
                    assembly_id = max_assembly_id + 1

                max_assembly_id = max(max_assembly_id, assembly_id)

                if "(" in assembly["oper_expression"] or "P" in assembly["oper_expression"]:
                    symmetry_ids = assembly["oper_expression"]
                else:
                    symmetry_ids = assembly["oper_expression"].split(",")

                chain_symmetry_seen = set()

                for pdb_chain in pdb_chains:
                    if pdb_chain in pdb_chain_to_author_chain:
                        author_chain = pdb_chain_to_author_chain[pdb_chain]
                        self.logger.info("Mapped pdb_chain %s to author_chain %s" % (pdb_chain,author_chain))
                        if not author_chain:
                            self.logger.info("author_chain is null, so skip this one")
                            continue
                    else:
                        self.logger.info("No mapping to author chain for pdb_chain %s" % pdb_chain)
                        continue

                    data = {}
                    data["pdb_id"] = pdb
                    data["chain_name"] = author_chain
                    data["assembly_id"] = assembly_id

                    if type(symmetry_ids) == list:
                        for symmetry_id in symmetry_ids:
                            chain_symmetry = author_chain + "&" + symmetry_id
                            if not chain_symmetry in chain_symmetry_seen:
                                # only add each author chain once, even if multiple different pdb chains for it
                                chain_symmetry_seen.add(chain_symmetry)
                                data["symmetry_id"] = symmetry_id
                                data["symmetry"] = symmetry_id_to_name[symmetry_id]

                                self.logger.info(data)

                                yield(mod.AssemblyInfo(**data))
                    else:
                        # complex operator expression, just store that as the name of the symmetry
                        data["symmetry_id"] = '0'
                        data["symmetry"] = symmetry_ids

                        self.logger.info(data)

                        yield(mod.AssemblyInfo(**data))
        except:
            try:
                ihm = cif_access.table('ihm_struct_assembly_details').rows
                for row in ihm:
                    self.logger.info('ihm row %s' % str(row))

                    assembly_id = row["assembly_id"]
                    # if not assembly_id in ["1",1]:
                    #     # seems to be working OK
                    #     self.logger.info('ihm has assembly other than 1, check it')

                    pdb_chain = row["asym_id"]
                    if not pdb_chain_to_author_chain:
                        author_chain = pdb_chain
                    elif pdb_chain in pdb_chain_to_author_chain:
                        author_chain = pdb_chain_to_author_chain[pdb_chain]
                        self.logger.info("Mapped pdb_chain %s to author_chain %s" % (pdb_chain,author_chain))
                    else:
                        self.logger.info("No mapping to author chain for pdb_chain %s" % pdb_chain)
                        continue

                    data = {}
                    data["pdb_id"] = pdb
                    data["assembly_id"] = assembly_id
                    data["chain_name"] = author_chain
                    data["symmetry_id"] = '1'
                    data["symmetry"] = '1_555'
                    yield(mod.AssemblyInfo(**data))
            except:
                self.logger.info('Not able to get assembly information for %s' % pdb)
                raise core.Skip("Not able to get assembly information, try again later")

