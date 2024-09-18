"""
Get PDB chain sequences that have not been mapped to Rfam.
Run cmsearch to map them.
More.

"""

from collections import defaultdict
import datetime
import gzip
import os

from pymotifs.constants import DATA_FILE_DIRECTORY
import pymotifs.core as core
from pymotifs import models as mod
from sqlalchemy import or_

# constants that should maybe be defined elsewhere
RFAM_ALIGNMENT_DIRECTORY = '/usr/local/pipeline/alignments'
INFERNAL_LOCATION = '/usr/local/pipeline/alignments/infernal-1.1.4'
INFERNAL_LOCATION = '/usr/local/pipeline/alignments/infernal-1.1.5'
MIN_PDB_CHAIN_LENGTH = 30

# lists of Rfam families to include in cross-domain alignments
# align according to the covariance model for the first family
joint_alignments = {}
joint_alignments['SSU'] = ['RF01960','RF01959','RF00177','RF02545','RF02542']
joint_alignments['LSU'] = ['RF02543','RF02540','RF02541','RF02546']

class Loader(core.Loader):
    merge_data = True
    mark = False
    allow_no_data = True

    # need to add stages that produce chain_info
    dependencies = set([])


    def read_mapping(self,filename):
        """
        Read a file like Rfam.pdb and store the chain to Rfam mappings there
        Note: only one mapping is stored per chain, and that is the one with the highest bit_score
              That should help with situations like 6YDP where one chain is too long and mapped to many families
        rfam_acc	pdb_id	chain	pdb_start	pdb_end	bit_score	evalue_score	cm_start	cm_end	hex_colour
        RF00001	4u3u	7	1	121	88.80	6.5e-20	1	119	8484c0
        """

        if not os.path.exists(filename):
            return {}

        if filename.endswith('.gz'):
            with gzip.open(filename,'rt') as f:
                lines = f.readlines()
        else:
            with open(filename,'r') as f:
                lines = f.readlines()

        chain_to_range_to_mapping = {}

        # work through entries by PDB id and chain to facilitate reading
        #for line in sorted(lines, key=lambda x: (x.split("\t")[1].upper(),x.split("\t")[2])):

        # work through entries one by one
        self.logger.info('Found %d lines in %s' % (len(lines),filename))
        count = 0
        for line in lines:
            if line.startswith('RF'):
                fields = line.replace('\n','').split("\t")
                if int(fields[3]) < int(fields[4]) or (fields[3] == '0' and fields[4] == '0'):
                    # only use positive strand matches, where positions are in order
                    chain = "%s|1|%s" % (fields[1].upper(),fields[2])  # like 4V9F|1|0

                    if not chain in chain_to_range_to_mapping:
                        chain_to_range_to_mapping[chain] = {}

                    header_list = ['rfam_acc','pdb_id','chain','pdb_start','pdb_end','bit_score','evalue_score','cm_start','cm_end','hex_colour']
                    md = {}  # mapping dictionary
                    for i,header in enumerate(header_list):
                        md[header] = fields[i]
                    md['line'] = line

                    range_string = '%09s-%09s' % (fields[3],fields[4])  # start and stop positions in the chain

                    chain_to_range_to_mapping[chain][range_string] = md

                    count += 1

                    """
                    # it would be possible to only keep one mapping for each chain
                    # that makes it easier to attach one name to a chain
                    # but for alignments, it's better to keep all mappings
                    if chain in chain_to_range_to_mapping:
                        if float(md['bit_score']) > float(chain_to_range_to_mapping[chain]['bit_score']):
                            # replace with mapping having higher score
                            self.logger.info('Replacing mapping of %-12s to %s with %s because it has higher bit score' % (chain,chain_to_range_to_mapping[chain]['rfam_acc'],md['rfam_acc']))
                            chain_to_range_to_mapping[chain] = md
                    else:
                        chain_to_range_to_mapping[chain] = md
                    """

        self.logger.info('Found %d mappings in %s' % (count,filename))

        return chain_to_range_to_mapping

    def write_mapping(self,filename,chain_to_range_to_mapping):
        """
        Write a file with a format like Rfam.pdb
        If hex_colour is not defined, substitute "rna3dhub" so it's clear where the mapping was done

        rfam_acc	pdb_id	chain	pdb_start	pdb_end	bit_score	evalue_score	cm_start	cm_end	hex_colour
        RF00001	4u3u	7	1	121	88.80	6.5e-20	1	119	8484c0
        """

        lines = []
        for chain in chain_to_range_to_mapping.keys():
            for range_string in chain_to_range_to_mapping[chain].keys():
                cm = chain_to_range_to_mapping[chain][range_string]

                if "line" in cm:
                    lines.append(cm['line'])
                else:
                    header_list = ['rfam_acc','pdb_id','chain','pdb_start','pdb_end','bit_score','evalue_score','cm_start','cm_end','hex_colour']
                    line = ''
                    for header in header_list:
                        if header in cm:
                            line += str(cm[header])
                        elif header == 'hex_coulour':
                            line += 'rna3dhub'

                        if header == 'hex_colour':
                            line += '\n'
                        else:
                            line += '\t'
                    lines.append(line)

        # write lines by Rfam family, increasing score
        lines = sorted(lines, key=lambda x: (x.split("\t")[0],float(x.split("\t")[5])))

        if filename.endswith('.gz'):
            with gzip.open(filename,'w') as f:
                f.writelines(lines)
        else:
            with open(filename,'w') as f:
                f.writelines(lines)

        return


    def write_pdb_chain_sequences_to_align(self,directory,rfam_family,chain_to_range_to_mapping,chain_to_sequence):
        """
        For this Rfam family, write the PDB chain sequences that are mapped to it.
        For a joint alignment, write a bunch of PDB chain sequences.
        Code will work on a single Rfam family or a joint alignment like 'LSU,RF02543,RF02540,RF02541,RF02546'
        """

        sequence_to_chain_range = {}

        for chain in chain_to_range_to_mapping.keys():
            for range_string in chain_to_range_to_mapping[chain].keys():
                cm = chain_to_range_to_mapping[chain][range_string]
                if cm['rfam_acc'] in rfam_family:
                    pdb_id = cm['pdb_id'].upper()
                    chain_name = cm['chain']
                    pdb_start = int(cm['pdb_start'])
                    pdb_end = int(cm['pdb_end'])

                    # previously-used chain_range notation
                    chain_range = "%s_%s_%s_%s" % (pdb_id,chain_name,pdb_start,pdb_end)

                    if chain in chain_to_sequence:
                        sequence = chain_to_sequence[chain]

                        if pdb_end <= len(sequence):
                            # extract only the part of the sequence that matches the model
                            sequence = sequence[pdb_start-1:pdb_end]

                            if not sequence in sequence_to_chain_range:
                                sequence_to_chain_range[sequence] = set()

                            sequence_to_chain_range[sequence].add(chain_range)

                    else:
                        self.logger.info('No sequence found for chain %s' % chain)
                        print('No sequence found for chain %s' % chain)

        num_chains = len(sequence_to_chain_range.keys())
        filename = os.path.join(directory,'sequences','%s_PDB_chains.fa' % rfam_family.split(",")[0])

        if num_chains > 0:
            with open(filename,'w') as f:
                for sequence in sequence_to_chain_range.keys():
                    f.write(">" + ",".join(sorted(sequence_to_chain_range[sequence])) + "\n")
                    f.write(sequence + "\n")
        else:
            # not sure why, but some empty files were written
            if os.path.exists(filename):
                os.remove(filename)
                self.logger.info('Removed empty file %s' % filename)
                print('Removed empty file %s' % filename)

        return num_chains


    def read_cmsearch_results(self,filename,rfam_family_to_minima):
        """
        Read the -tblout file from cmsearch
        Save each line, then keep non-overlapping results for each chain
        """

        if not os.path.exists(filename):
            return {}

        if filename.endswith('.gz'):
            with gzip.open(filename,'r') as f:
                lines = f.readlines()
        else:
            with open(filename,'r') as f:
                lines = f.readlines()

        results = []
        chain_to_range_to_mapping = {}
        chain_to_positions_taken = {}
        for line in lines:
            if not line.startswith("#"):
                fields = line.replace('\n','').split()
                if not len(fields) == 18:
                    print('Unexpected number of fields in cmsearch results line: %s' % line)
                chains = fields[0]
                for chain in chains.split(","):
                    # prepare for each chain, leave a mark for chains not mapped
                    chain_to_range_to_mapping[chain] = {}
                    chain_to_positions_taken[chain] = set()

                    result = {}
                    result['chain'] = chain
                    result['query_name'] = fields[2]
                    result['accession'] = fields[3]
                    result['mdl'] = fields[4]
                    result['mdl_from'] = fields[5]
                    result['mdl_to'] = fields[6]
                    result['seq_from'] = fields[7]
                    result['seq_to'] = fields[8]
                    result['strand'] = fields[9]
                    result['trunc'] = fields[10]
                    result['pass'] = fields[11]
                    result['gc'] = fields[12]
                    result['bias'] = fields[13]
                    result['score'] = float(fields[14])
                    result['E-value'] = fields[15]
                    result['inc'] = fields[16]
                    result['description'] = ' '.join(fields[17:])

                    # only use positive strand matches
                    if result['strand'] == '+':
                        # check other criteria for a good mapping
                        length = float(result['seq_to'])-float(result['seq_from'])
                        score_per_length = result['score'] / (0.01+length)  # score per length of match
                        rfam_family = result['accession']

                        result['length'] = length
                        result['score_per_length'] = score_per_length

                        if rfam_family in rfam_family_to_minima:
                            # have a higher score than Rfam maps, and have a decent score per length of match
                            if result['score'] >= rfam_family_to_minima[rfam_family]['score'] and score_per_length >= min(0.8,rfam_family_to_minima[rfam_family]['score_per_length']):
                                # new result is better than existing result
                                results.append(result)
                        else:
                            if result['score'] > 30 and score_per_length > 0.3 and length > MIN_PDB_CHAIN_LENGTH:
                                results.append(result)
                            elif score_per_length > 0.8:
                                results.append(result)

        # go from best match to worst match, keep non-overlapping matches for each chain
        results = sorted(results, key = lambda x: (x['score']), reverse=True)
        for result in results:
            chain = result['chain']
            new_range = set(range(int(result['seq_from']),int(result['seq_to'])))
            if not new_range & chain_to_positions_taken[chain]:
                # no overlap with previous mappings
                chain_to_positions_taken[chain] |= new_range
                range_string = '%09s-%09s' % (result['seq_from'],result['seq_to'])  # start and stop positions in the chain
                chain_to_range_to_mapping[chain][range_string] = result

                if len(chain_to_range_to_mapping[chain]) > 1:
                    for range_string in sorted(chain_to_range_to_mapping[chain].keys(), key = lambda x: (chain_to_range_to_mapping[chain][x]['score']), reverse=True):
                        result = chain_to_range_to_mapping[chain][range_string]
                        self.logger.info('Mapping %-12s positions %s to %s score %7.2f' % (chain,range_string,result['accession'],result['score']))
                        print('Mapping %-12s positions %s to %s score %7.2f' % (chain,range_string,result['accession'],result['score']))

        return chain_to_range_to_mapping


    def to_process(self, pdbs, **kwargs):
        """
        This is the starting point for map_to_rfam.
        Ignore pdbs because we want to map obsolete pdbs as well for backward compatibility.
        """

        print("Starting to_process for map_to_rfam")

        # track rfam families with new PDB chains that need to be aligned
        rfam_families_to_align = set()

        # read the current list of Rfam families that need to be aligned
        # if the previous process failed, it will try again now
        rfam_to_align_filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'rfam_families_to_align.txt')
        if os.path.exists(rfam_to_align_filename):
            with open(rfam_to_align_filename,'r') as f:
                for line in f.readlines():
                    rfam_families_to_align.add(line.replace('\n',''))

        # read Rfam.pdb.gz to identify chains that Rfam has mapped; Rfam.pdb.gz needs to be downloaded occasionally
        self.logger.info('Reading Rfam mappings')
        chain_to_range_to_mapping = self.read_mapping(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'Rfam.pdb.gz'))

        # find minimum score in each Rfam family according to Rfam mapping
        self.logger.info('Getting min Rfam scores in each family')
        rfam_family_to_minima = {}
        for chain, range_to_mapping in chain_to_range_to_mapping.items():
            for range_string, rfam_md in range_to_mapping.items():
                rfam_family = rfam_md['rfam_acc']
                score = float(rfam_md['bit_score'])
                length = float(rfam_md['pdb_end'])-float(rfam_md['pdb_start'])
                score_per_length = score / (0.01+length)  # score per length of match
                if not rfam_family in rfam_family_to_minima:
                    rfam_family_to_minima[rfam_family] = {'score':score,'length':length,'score_per_length':score_per_length}
                else:
                    if score < rfam_family_to_minima[rfam_family]['score']:
                        rfam_family_to_minima[rfam_family]['score'] = score
                    if length < rfam_family_to_minima[rfam_family]['length']:
                        rfam_family_to_minima[rfam_family]['length'] = length
                    if score_per_length < rfam_family_to_minima[rfam_family]['score_per_length']:
                        rfam_family_to_minima[rfam_family]['score_per_length'] = score_per_length

        # read pdb_chain_to_rfam.txt to identify chains that are already mapped by Rfam or rna3dhub
        self.logger.info('Reading previous mappings')
        rna3dhub_chain_to_range_to_mapping = self.read_mapping(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'pdb_chain_to_rfam.txt'))

        self.logger.info('rna3dhub_chain_to_range_to_mapping has %i chains' % len(rna3dhub_chain_to_range_to_mapping.keys()))

        # if a chain is not mapped by Rfam but was mapped here, add it
        self.logger.info('chain_to_range_to_mapping has %i chains' % len(chain_to_range_to_mapping.keys()))
        self.logger.info('Adding chains from pdb_chain_to_rfam.txt')
        for chain in rna3dhub_chain_to_range_to_mapping.keys():
            if not chain in chain_to_range_to_mapping:
                chain_to_range_to_mapping[chain] = rna3dhub_chain_to_range_to_mapping[chain]
        self.logger.info('chain_to_range_to_mapping has %i chains' % len(chain_to_range_to_mapping.keys()))

        # retrieve from the database all chains and sequences that can reasonably map to Rfam
        # or just retrieve all of them; the short ones don't amount to much
        # This avoids using pdbseqres
        self.logger.info('Getting RNA chains from the database')
        with self.session() as session:
            query = session.query(mod.ChainInfo.pdb_id,
                                  mod.ChainInfo.chain_name,
                                  mod.ChainInfo.sequence,
                                  mod.ChainInfo.chain_length,
                                  mod.PdbInfo.title).\
                outerjoin(mod.PdbInfo, mod.ChainInfo.pdb_id == mod.PdbInfo.pdb_id).\
                filter(or_(mod.ChainInfo.entity_macromolecule_type.contains('Polyribonucleotide'),
                           mod.ChainInfo.entity_macromolecule_type.contains('Hybrid')))
                #filter(mod.ChainInfo.entity_macromolecule_type.contains('Polyribonucleotide'))
                #filter(mod.ChainInfo.chain_length > 1000).\

        chain_to_sequence = {}
        sequence_to_chains = defaultdict(list)  # many chains could have the same sequence
        chain_to_title = {}

        for r in query:
            chain = "%s|1|%s" % (r.pdb_id.upper(),r.chain_name)
            chain_to_sequence[chain] = r.sequence
            chain_to_title[chain] = r.title

            if len(r.sequence) > MIN_PDB_CHAIN_LENGTH:
                sequence_to_chains[r.sequence].append(chain)

        # if any unmapped sequences are the same as a chain that is already mapped, map them like before
        # this avoids having to run cmsearch on them, and promotes consistency.  Also propagates any errors!
        self.logger.info('Identifying chains whose sequences are already mapped to Rfam')
        for sequence, chains in sequence_to_chains.items():
            new_chains = set(chains) - set(chain_to_range_to_mapping.keys())
            mapped_chains = set(chains) & set(chain_to_range_to_mapping.keys())

            if len(mapped_chains) > 1:
                # it's conceivable that the same sequence has different mapping information
                # do some simple checks to see if they are the same
                range_to_count = defaultdict(int)
                for chain in mapped_chains:
                    for range_string in chain_to_range_to_mapping[chain].keys():
                        range_to_count[range_string] += 1
                counts = list(range_to_count.values())
                if min(counts) < max(counts):
                    self.logger.info('Chains %s have different mappings' % sorted(mapped_chains))
                    print('Chains %s have different mappings' % sorted(mapped_chains))

            # if there are both new chains and already mapped chains with this sequence
            if len(new_chains) > 0 and len(mapped_chains) > 0:
                # just pick one to work with
                mapped_chain = list(mapped_chains)[0]
                for chain in new_chains:

                    if not chain in chain_to_range_to_mapping:
                        chain_to_range_to_mapping[chain] = {}
                    for range_string in chain_to_range_to_mapping[mapped_chain].keys():
                        cm = {}
                        for key, value in chain_to_range_to_mapping[mapped_chain][range_string].items():
                            if not key == 'line':
                                cm[key] = value
                        cm['pdb_id'] = chain.split('|')[0]
                        cm['chain'] = chain.split('|')[2]
                        chain_to_range_to_mapping[chain][range_string] = cm

                    rfam_mapped_to = [chain_to_range_to_mapping[chain][range_string]['rfam_acc'] for range_string in chain_to_range_to_mapping[chain].keys()]
                    print('Chain %s has same sequence as %s which maps to %s' % (chain,mapped_chain,rfam_mapped_to))
                    self.logger.info('Chain %s has same sequence as %s which maps to %s' % (chain,mapped_chain,rfam_mapped_to))
                    # self.logger.info('New dictionary is %s' % chain_to_range_to_mapping[chain])

                    # PDB chains for this family need to be aligned, new combined alignment produced
                    rfam_families_to_align |= set(rfam_mapped_to)

        self.logger.info('chain_to_range_to_mapping has %i chains' % len(chain_to_range_to_mapping.keys()))

        # accumulate unmapped chain sequences, write to a fasta file, run cmsearch
        lines = []
        chains_to_map = set()
        # print(sorted(chain_to_range_to_mapping.keys()))

        for sequence, chains in sequence_to_chains.items():
            #self.logger.info('Looking for %s in the mapped chains' % chains[0])
            #print('Looking for %s in the mapped chains' % chains[0])
            if not chains[0] in chain_to_range_to_mapping:
                chains_to_map |= set(chains)
                lines.append('>' + ','.join(chains))
                lines.append(sequence)
                self.logger.info('Did not find %s in the mapped chains' % chains[0])
                print('Did not find %s in the mapped chains' % chains[0])

        if len(lines) > 0:
            print('Found %i unmapped sequences, writing to chains_for_cmsearch.fa' % (len(lines)/2))
            self.logger.info('Found %i unmapped sequences, writing to chains_for_cmsearch.fa' % (len(lines)/2))
            with open(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'chains_for_cmsearch.fa'),'w') as f:
                f.write('\n'.join(lines))

            # run cmsearch on chains_to_map, save results in cmsearch_results.txt
            print('Running cmsearch on chains_for_cmsearch.fa')
            self.logger.info('Running cmsearch on chains_for_cmsearch.fa')
            cmsearch_location = os.path.join(INFERNAL_LOCATION,'src','cmsearch')
            fasta_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'chains_for_cmsearch.fa')
            output_table = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'cmsearch_results.txt')
            rfam_cm_models = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'Rfam.cm.gz')
            os.system('%s --cpu 4 --tblout %s %s %s' % (cmsearch_location,output_table,rfam_cm_models,fasta_file))

        # read results of cmsearch, evaluate to be good enough, take non-overlapping mappings
        chain_to_range_to_cmsearch_results = self.read_cmsearch_results(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'cmsearch_results.txt'),rfam_family_to_minima)

        # accumulate new data to write to pdb_chain_to_rfam.txt
        self.logger.info("Writing any new lines to pdb_chain_to_rfam.txt")
        output_lines = []
        for chain in chains_to_map:
            if chain in chain_to_range_to_mapping:
                # already mapped by Rfam or in a previous week, so skip
                continue
            elif chain in chain_to_range_to_cmsearch_results and len(chain_to_range_to_cmsearch_results[chain]) > 0:
                # loop over ranges for that chain
                range_to_result = chain_to_range_to_cmsearch_results[chain]

                for range_string in range_to_result.keys():
                    result = range_to_result[range_string]

                    if not chain in chain_to_range_to_mapping:
                        chain_to_range_to_mapping[chain] = {}

                    output_line = '%s %7.2f %7.2f %4s %4s %-12s %s\n' % (result['accession'],result['score'],result['score_per_length'],result['seq_from'],result['seq_to'],chain,chain_to_title[chain])
                    output_lines.append(output_line)

                    # re-format the data for pdb_chain_to_rfam.txt
                    # ['rfam_acc','pdb_id','chain','pdb_start','pdb_end','bit_score','evalue_score','cm_start','cm_end','hex_colour']
                    new_md = {}
                    new_md['rfam_acc'] = result['accession']
                    fields = chain.split('|')

                    new_md['pdb_id'] = fields[0]
                    new_md['chain'] = fields[2]
                    new_md['pdb_start'] = result['seq_from']
                    new_md['pdb_end'] = result['seq_to']
                    new_md['bit_score'] = result['score']
                    new_md['evalue_score'] = result['E-value']
                    new_md['cm_start'] = result['mdl_from']
                    new_md['cm_end'] = result['mdl_to']
                    new_md['hex_colour'] = 'rna3dhub'

                    range_string = '%09s-%09s' % (result['seq_from'],result['seq_to'])  # start and stop positions in the chain
                    chain_to_range_to_mapping[chain][range_string] = new_md

                    # record that the Rfam family needs to be re-aligned
                    rfam_families_to_align.add(result['accession'])

                    # if the Rfam family is part of a joint alignment (e.g., SSU or LSU cross domain),
                    # add the whole group of families like 'LSU,RF02543,RF02540,RF02541,RF02546'
                    for molecule in joint_alignments.keys():
                        if result['accession'] in joint_alignments[molecule]:
                            rfam_families_to_align.add("%s,%s" % (molecule,",".join(joint_alignments[molecule])))

                    #print('Mapping %-12s to %s' % (chain,result['accession']))
                    self.logger.info('Mapping %-12s to %s' % (chain,result['accession']))

            else:
                # record that the chain was considered but not mapped
                # this includes chains which cmsearch simply skips for whatever reason, like having N in the sequence
                new_md = {}
                new_md['rfam_acc'] = 'RF00000'
                fields = chain.split('|')
                new_md['pdb_id'] = fields[0]
                new_md['chain'] = fields[2]
                new_md['pdb_start'] = 0
                new_md['pdb_end'] = 0
                new_md['bit_score'] = 0.0
                new_md['evalue_score'] = 0.0
                new_md['cm_start'] = 0
                new_md['cm_end'] = 0
                new_md['hex_colour'] = 'rna3dhub'

                if not chain in chain_to_range_to_mapping:
                    chain_to_range_to_mapping[chain] = {}
                range_string = '%09s-%09s' % (0,0)  # pretend start and stop positions in the chain
                chain_to_range_to_mapping[chain][range_string] = new_md

                print('No mapping found for %-12s' % (chain))
                self.logger.info('No mapping found for %-12s' % (chain))

        # write information about the decisions
        current_datetime = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        with open(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'mapping_decision','mapping_decision_%s.txt' % current_datetime),'w') as f:
            # write basic information about each Rfam family's mappings
            for rfam_family, min_score in sorted(rfam_family_to_minima.items()):
                f.write('%s Min score %7.2f Min length %5d Min ratio %7.2f\n' % (rfam_family,min_score['score'],min_score['length'],min_score['score_per_length']))
            for output_line in sorted(output_lines):
                f.write(output_line)

        # write all data to pdb_chain_to_rfam.txt
        print('Writing mappings to pdb_chain_to_rfam.txt')
        self.logger.info('chain_to_range_to_mapping has %i chains' % len(chain_to_range_to_mapping.keys()))
        self.write_mapping(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'pdb_chain_to_rfam.txt'),chain_to_range_to_mapping)

        # when you need to re-align everything, you can use this list
        if False:
            rfam_families_to_align = set()
            for chain in chain_to_range_to_mapping.keys():
                for range_string in chain_to_range_to_mapping[chain].keys():
                    rfam_families_to_align.add(chain_to_range_to_mapping[chain][range_string]['rfam_acc'])

        # work just with specific Rfam families
        # rfam_families_to_align = set(['RF00023'])  # very slow to align all the sequences!
        # rfam_families_to_align = set(['RF01959'])  # good to practice with long alignment
        #rfam_families_to_align = set(['RF00005','RF02543'])

        # practice with SSU and LSU alignments
        # rfam_families_to_align = set()
        # for molecule in joint_alignments.keys():
        #     rfam_families_to_align.add("%s,%s" % (molecule,",".join(joint_alignments[molecule])))

        sorted_rfam_families_to_align = sorted(rfam_families_to_align)

        # write a list of Rfam families that have new chains and so need to be re-aligned
        print('Writing rfam_families_to_align.txt')
        self.logger.info('Writing rfam_families_to_align.txt')
        # only ever append to this list here
        # the only thing that can take away from the list is the aligning program
        with open(rfam_to_align_filename,'w') as f:
            for rfam_family in sorted_rfam_families_to_align:
                if not 'RF00000' in rfam_family:
                    # write PDB chain sequences to align in this Rfam family
                    num_chains = self.write_pdb_chain_sequences_to_align(RFAM_ALIGNMENT_DIRECTORY,rfam_family,chain_to_range_to_mapping,chain_to_sequence)
                    print('Wrote %4d PDB chain sequences to align to %s' % (num_chains,rfam_family))
                    self.logger.info('Wrote %4d PDB chain sequences to align to %s' % (num_chains,rfam_family))
                    f.write('%s\n' % rfam_family)

        if len(sorted_rfam_families_to_align) == 0:
            print('No new Rfam families to align')
            self.logger.info('No new Rfam families to align')
            sorted_rfam_families_to_align = ['RF00000']

        return sorted_rfam_families_to_align

    def has_data(self, pdb_dict, **kwargs):
        """
        This method can query the database after to_process but before data method
        to see if there is already data in the database, and if so, it returns True,
        and the data method does not have to work on this item.
        """

        return False

    def remove(self, pdb, **kwargs):

        return True


    def extract_covariance_models(self,rfam_families):
        """
        Extract Infernal and HMMER models of Rfam families with a 3D structure
        """

        rfam_model_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'Rfam.cm.gz')

        with gzip.open(rfam_model_file,'r') as file:
            cm_data = file.read().split("//")

        for i in range(len(cm_data)):
            model = cm_data[i].strip()
            if model.startswith("INFERNAL"):
                for family in rfam_families:
                    if family in model:
                        family_model_filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'covariance_models','%s.cm' % family)
                        with open(family_model_filename, 'w') as out_file:
                            out_file.write(model+"\n//"+cm_data[i+1]+"//")
                            print("Extracted %s.cm" % family)


    def rectify_sequence(self,sequence):
        """
        Change . to -
        Change leading and trailing - to .
        """
        sequence = sequence.replace(".","-")

        j = 0
        while j < len(sequence) and sequence[j] == '-':
            j += 1

        k = len(sequence)-1
        while k >= 0 and sequence[k] == '-':
            k -= 1

        new_sequence = '.'*j + sequence[j:k+1] + '.'*(len(sequence)-k-1)

        if len(new_sequence) == len(sequence):
            pass
        else:
            print('Old length %i new length %i' % (len(sequence),len(new_sequence)))

        return new_sequence

    def rectify_alignment_file(self,alignment_file,alignment_file_out_gz=None):
        """
        Read alignment in various formats, write out in fasta format
        cmalign sto format splits long sequences onto multiple lines; combine
        cmalign AFA format splits long sequences onto multiple lines; combine
        cmalign Pfam format has long header on one line, then short header and sequence later
        Also, replace leading and trailing - with .
        """

        if not alignment_file_out_gz:
            alignment_file_out_gz = alignment_file+".gz"

        new_lines = []
        if ".fa" in alignment_file.lower():
            # fasta format, not so bad to read all of the lines at once
            with open(alignment_file,'r') as f:
                lines = f.readlines()

            i = 0
            while i < len(lines):
                if lines[i].startswith(">"):
                    new_lines.append(lines[i])
                    i += 1
                    sequence = ""
                    while i < len(lines) and not lines[i].startswith(">"):
                        sequence += lines[i].replace('\n','')
                        i += 1

                    new_sequence = self.rectify_sequence(sequence)
                    new_lines.append(new_sequence + "\n")

        elif ".sto" in alignment_file.lower():
            # Infernal Stockholm format for this problem
            id_to_header = {}
            id_to_sequence = {}
            id_list_sequence = []
            id_list_pdb = []

            i = 0
            with open(alignment_file,'r') as f:
                for line in f:
                    # remove trailing \n
                    line = line.rstrip('\n')
                    if line.startswith("#=GS"):
                        # map first field of name to rest of name
                        fields = line.split()
                        # fields[0] is "#=GS" which we don't need
                        # fields[1] is typically the accession number
                        id = fields[1]
                        # fields[2] is "DE" which we don't need
                        # fields[3] is taxid in our sequence data
                        # fields[4] is domain in our sequence data
                        # fields[5:] is species and description
                        header = id + " " + " ".join(fields[3:])  # put space between accession and taxid
                        id_to_header[id] = header
                        id_to_sequence[id] = ""
                        id_list_sequence.append(id)           # keep lines in order
                    elif line.startswith("#") or line.startswith("/"):
                        pass
                    elif len(line) > 0:
                        fields = line.split()
                        if len(fields) == 2:
                            id = fields[0]
                            if not id in id_to_header:
                                id_to_header[id] = id    # PDB chains are like this
                                id_to_sequence[id] = ""  # first time being seen
                                id_list_pdb.append(id)      # keep lines in order
                            id_to_sequence[id] += fields[1]
                        else:
                            print("Line %d which is %s has %i fields" % (i,line,len(fields)))

                    i += 1

            for id in id_list_pdb + id_list_sequence:
                new_lines.append(">%s\n" % id_to_header[id])
                sequence = id_to_sequence[id]
                new_sequence = self.rectify_sequence(sequence)
                new_lines.append(new_sequence + "\n")

        elif ".pfam" in alignment_file.lower():
            # Pfam format, from how I have heard it described, pretty much Stockholm
            k = 0
            with open(alignment_file,'r') as f:
                for line in f:
                    if line.startswith("#") or line.startswith("/") or len(line) < 5:
                        continue
                    if k == 0:
                        # find the first column of sequence, using the first line of the file
                        k = len(line)-1
                        while k > 0 and line[k] != ' ':
                            k -= 1
                    new_lines.append(">"+line[:k].rstrip() + "\n")  # remove spaces after header
                    new_sequence = self.rectify_sequence(line[k+1:].replace('\n',''))
                    new_lines.append(new_sequence + "\n")

        with gzip.open(alignment_file_out_gz,'wt') as f:
            f.writelines(new_lines)


    def align_large_sequence_family(self,rfam_family):
        """
        Align sequences in this family progressively, a few thousand at a time
        """

        cmalign_location = os.path.join(INFERNAL_LOCATION,'src','cmalign')
        alimerge_location = os.path.join(INFERNAL_LOCATION,'easel','miniapps','esl-alimerge')
        model_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'covariance_models','%s.cm' % rfam_family)

        # read the entire sequence file
        sequence_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_rfam_sequences.fa.gz' % rfam_family)
        with gzip.open(sequence_file,'r') as f:
            lines = f.readlines()

        self.logger.info('Found %i lines in large sequence family %s' % (len(lines),sequence_file))
        print('Found %i lines in large sequence family %s' % (len(lines),sequence_file))

        start_index = 0

        if rfam_family == 'RF00005':
            lines_delta = 100000
        elif rfam_family == 'RF02543':
            lines_delta = 200
        else:
            lines_delta = int(len(lines)/20.0)

        # small temporary file for alignment of a portion of the sequences
        sequence_temp = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_rfam_sequences_temp.fa' % rfam_family)

        # growing alignment file
        alignment_growing = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_rfam_sequences_growing.sto' % rfam_family)

        # final alignment file
        seq_alignment_file_sto = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_rfam_sequences.sto' % rfam_family)

        while start_index < len(lines):
            last_index = min(start_index+lines_delta,len(lines))
            if len(lines) - last_index < lines_delta*0.4:
                # just do all the rest of the sequences
                last_index = len(lines)
                lines_delta = last_index - start_index

            with open(sequence_temp,'w') as f:
                f.writelines(lines[start_index:last_index])

            # first time write to the growing alignment file, after that write to a temp file
            if start_index == 0:
                alignment_temp = alignment_growing
            else:
                alignment_temp = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_rfam_sequences_temp.sto' % rfam_family)

            # run infernal; default mxsize is 128, 1200 covers most big alignments
            # mxsize 1000 worked better for LSU_PDB_chains on 2024-07-08
            command = '%s --mxsize 1000 --noprob --outformat Pfam %s %s > %s' % (cmalign_location,model_file,sequence_temp,alignment_temp)
            self.logger.info('Start %5d last %5d running command: %s' % (start_index,last_index,command))
            print('Start %5d last %5d running command: %s' % (start_index,last_index,command))
            os.system(command)

            # check to see that Infernal really produced an alignment
            # if not, remove and bail out
            if os.path.exists(alignment_temp) and os.path.getsize(alignment_temp) == 0:
                os.remove(alignment_temp)
                print('  Infernal failed to align %s' % sequence_temp)
                self.logger.info('Command was: %s' % command)
                raise core.Skip("Something went wrong with aligning sequences for %s; try again next time" % (rfam_family))

            # merge the temporary new alignment with the growing full alignment
            if start_index > 0:
                # move the previously produced large alignment so it can be added to
                command = 'mv %s %s' % (seq_alignment_file_sto,alignment_growing)
                os.system(command)

                command = '%s --outformat pfam -o %s %s %s' % (alimerge_location,seq_alignment_file_sto,alignment_temp,alignment_growing)
                print('Combining alignments')
                print('Running command: %s' % command)
                os.system(command)

            start_index += lines_delta


        os.remove(alignment_temp)
        os.remove(alignment_growing)



    def data(self, rfam_family, **kwargs):
        """
        Input is a list of Rfam families that need to be aligned.
        First we align the PDB chains.
        That way, the PDB chains get aligned quickly, and then we work on full families, which is slower.
        rfam_family could be a single family or a list for a joint alignment like 'LSU,RF02543,RF02540,RF02541,RF02546'
        """

        # read the current list of Rfam families that need to be aligned
        rfam_families_to_align = set()
        with open(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'rfam_families_to_align.txt'),'r') as f:
            for line in f.readlines():
                rfam_families_to_align.add(line.replace('\n',''))

        family_id = rfam_family.split(',')[0]  # e.g., LSU, RF02543

        if family_id == 'LSU':
            model_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'covariance_models','RF02543.cm') # eukaryotic LSU
        elif family_id == 'SSU':
            model_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'covariance_models','RF01960.cm') # eukaryotic SSU
        else:
            model_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'covariance_models','%s.cm' % rfam_family)

        cmalign_location = os.path.join(INFERNAL_LOCATION,'src','cmalign')
        alimerge_location = os.path.join(INFERNAL_LOCATION,'easel','miniapps','esl-alimerge')

        if not os.path.exists(model_file):
            # extract the covariance models that will be needed
            self.extract_covariance_models([rfam_family])

            if not os.path.exists(model_file):
                raise core.Skip("Covariance model for %s is not found" % (rfam_family))

        pdb_chain_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_PDB_chains.fa' % family_id)
        if not os.path.exists(pdb_chain_file):
            raise core.Skip("PDB chain file for %s is not found" % (rfam_family))

        # align PDB chains mapped to this Rfam family
        self.logger.info("Aligning PDB chains mapped to %s using Infernal cmalign" % rfam_family)
        print("Aligning PDB chains mapped to %s using Infernal cmalign" % rfam_family)

        # run infernal cmalign to align PDB chain sequences to the model
        # default mxsize is 128, 1200 covers most big alignments
        # LSU pdb chains got to 13.2 GB on 2024-09-12, about 10 minutes into a 13-minute run
        pdb_alignment_file_sto = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_PDB_chains.sto' % family_id)
        command = '%s --mxsize 1200 --noprob --outformat Pfam %s %s > %s' % (cmalign_location,model_file,pdb_chain_file,pdb_alignment_file_sto)
        self.logger.info('Running command: %s' % command)
        print('Running command: %s' % command)
        os.system(command)

        # even if Infernal fails, the .sto alignment file will be created, so then a .gz version will be made
        # check to see that Infernal really produced an alignment
        # if so gzip, if not remove
        if not os.path.exists(pdb_alignment_file_sto):
            self.logger.info('Alignment file %s does not exist' % (pdb_alignment_file_sto))
        elif os.path.getsize(pdb_alignment_file_sto) == 0:
            os.remove(pdb_alignment_file_sto)
            self.logger.info('Infernal failed to align %s' % pdb_chain_file)
            self.logger.info('Command was: %s' % command)
            raise core.Skip("Something went wrong with aligning PDB chain sequences for %s" % (rfam_family))
        else:
            # read PDB chain alignment file, put sequence all on one line, replace leading and trailing - with .
            alignment_file_fa = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_PDB_chains.fa' % family_id)
            self.rectify_alignment_file(pdb_alignment_file_sto,alignment_file_fa+".gz")

            # check to see if sequences from Rfam full family are already collected; if not, note that
            sequence_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_rfam_sequences.fa' % rfam_family)
            sequence_file_gz = sequence_file + ".gz"

            if not os.path.exists(sequence_file) and not os.path.exists(sequence_file_gz):
                # note that the sequence file is needed for this Rfam family
                rfam_sequence_file_needed = set()
                sequence_needed_filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'rfam_sequence_file_needed.txt')
                # read list of Rfam families that need sequence files
                if os.path.exists(sequence_needed_filename):
                    with open(sequence_needed_filename,'r') as f:
                        for line in f.readlines():
                            rfam_sequence_file_needed.add(line.replace('\n',''))
                rfam_sequence_file_needed.add(rfam_family)
                # write list of Rfam families that need sequence files
                with open(sequence_needed_filename,'w') as f:
                    for rf in sorted(rfam_sequence_file_needed):
                        f.write('%s\n' % rf)

                self.logger.info('No sequence file %s found for %s' % (sequence_file,rfam_family))
                print('No sequence file found for %s' % rfam_family)
                raise core.Skip("Rfam sequence file for %s is not found; check again on the next run" % (rfam_family))

            # check to see if sequences from Rfam full family are already aligned; if not, align them
            seq_alignment_file_sto = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_rfam_sequences.sto' % rfam_family)
            seq_alignment_file_sto_gz = seq_alignment_file_sto + ".gz"

            if not os.path.exists(seq_alignment_file_sto) and not os.path.exists(seq_alignment_file_sto_gz):
                self.logger.info("Aligning sequences  from %s full family using Infernal cmalign" % rfam_family)
                print("Aligning sequences  from %s full family using Infernal cmalign" % rfam_family)

                # some Rfam family have so many sequences, or such long sequences, that they need a different approach
                if rfam_family in ['RF00005','RF02543']:
                    self.align_large_sequence_family(rfam_family)

                else:
                    if not os.path.exists(sequence_file):
                        os.system("gunzip -f %s" % (sequence_file_gz))

                    # run infernal; default mxsize is 128, 1200 covers most big alignments
                    # command = '%s --mxsize 1200 --noprob %s %s > %s' % (cmalign_location,model_file,sequence_file,alignment_file)
                    command = '%s --mxsize 1200 --noprob --outformat Pfam %s %s > %s' % (cmalign_location,model_file,sequence_file,seq_alignment_file_sto)
                    print('Running command: %s' % command)
                    os.system(command)

                    # zip the raw sequence file again
                    os.system("gzip -f %s" % (sequence_file))

                    # even if Infernal fails, the alignment file will be created
                    # check to see that Infernal really produced an alignment
                    # if not, remove
                    if os.path.exists(seq_alignment_file_sto) and os.path.getsize(seq_alignment_file_sto) == 0:
                        os.remove(seq_alignment_file_sto)
                        self.logger.info('Infernal failed to align %s' % sequence_file)
                        print('Infernal failed to align %s' % sequence_file)
                        raise core.Skip("Something went wrong with aligning sequences for %s; try again next time" % (rfam_family))

            # merge the PDB chain alignment and the Rfam sequence alignment
            combined_sto = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_combined.sto' % family_id)

            if not os.path.exists(seq_alignment_file_sto):
                os.system("gunzip -f %s" % (seq_alignment_file_sto_gz))

            command = '%s --outformat pfam -o %s %s %s' % (alimerge_location,combined_sto,pdb_alignment_file_sto,seq_alignment_file_sto)
            self.logger.info('Combining PDB chain and Rfam sequence alignments')
            print('Combining PDB chain and Rfam sequence alignments')
            print('Running command: %s' % command)
            os.system(command)

            # read combined alignment file, put sequence all on one line, replace leading and trailing - with .
            self.logger.info("Rectifying alignment")
            print("Rectifying alignment")
            alignment_file_fa = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_combined.fa' % family_id)
            self.rectify_alignment_file(combined_sto,alignment_file_fa+".gz")

            # remove rfam_family from the list of rfam families to align
            with open(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'rfam_families_to_align.txt'),'w') as f:
                for rf in sorted(rfam_families_to_align):
                    # remove rfam_family from rfam_families_to_align
                    if not rf == rfam_family:
                        f.write('%s\n' % rf)

            # clean up files we don't need
            os.remove(combined_sto)
            os.remove(pdb_alignment_file_sto)
            os.system("gzip -f %s" % (seq_alignment_file_sto))

