"""
Get PDB chain sequences that have not been mapped to Rfam.
Run cmsearch to map them.
More.
"""

# imports that we really need
from collections import defaultdict
import datetime
import gzip
import os
import pymotifs.core as core
from pymotifs import models as mod
from sqlalchemy import or_

# constants that should maybe be defined elsewhere
RFAM_ALIGNMENT_DIRECTORY = '/usr/local/pipeline/alignments'
INFERNAL_LOCATION = '/usr/local/pipeline/alignments/infernal-1.1.5'
INFERNAL_LOCATION = '/usr/local/pipeline/alignments/infernal-1.1.4'
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
            with gzip.open(filename,'r') as f:
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

        print('write_pdb 1',rfam_family)

        for chain in chain_to_range_to_mapping.keys():
            for range_string in chain_to_range_to_mapping[chain].keys():
                cm = chain_to_range_to_mapping[chain][range_string]
                if cm['rfam_acc'] in rfam_family:
                    pdb_id = cm['pdb_id'].upper()
                    chain_name = cm['chain']
                    pdb_start = int(cm['pdb_start'])
                    pdb_end = int(cm['pdb_end'])

                    # print('write_pdb 2',rfam_family)


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
        print('Getting min Rfam criteria')
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

        self.logger.info('chain_to_range_to_mapping has %i chains' % len(chain_to_range_to_mapping.keys()))


        # accumulate unmapped chain sequences, write to a fasta file, run cmsearch
        lines = []
        # print(sorted(chain_to_range_to_mapping.keys()))

        for sequence, chains in sequence_to_chains.items():
            #self.logger.info('Looking for %s in the mapped chains' % chains[0])
            #print('Looking for %s in the mapped chains' % chains[0])
            if not chains[0] in chain_to_range_to_mapping:
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
            if False:
                cmsearch_location = os.path.join(INFERNAL_LOCATION,'src','cmsearch')
                fasta_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'chains_for_cmsearch.fa')
                output_table = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'cmsearch_results.txt')
                rfam_cm_models = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'Rfam.cm.gz')
                os.system('%s --cpu 4 --tblout %s %s %s' % (cmsearch_location,output_table,rfam_cm_models,fasta_file))

        # read results of cmsearch, evaluate to be good enough, take non-overlapping mappings
        chain_to_range_to_cmsearch_results = self.read_cmsearch_results(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'cmsearch_results.txt'),rfam_family_to_minima)

        # accumulate new data to write to pdb_chain_to_rfam.txt
        output_lines = []
        for chain in chain_to_range_to_cmsearch_results.keys():
            range_to_result = chain_to_range_to_cmsearch_results[chain]

            if chain in chain_to_range_to_mapping:
                # already mapped by Rfam or in a previous week, so skip
                continue
            elif len(range_to_result.keys()) > 0:
                # loop over ranges for that chain
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
                    rfam_families_to_align.add(result['accession']+"_combined")  # also need to align with sequences

                    # if the Rfam family is part of a joint alignment (e.g., SSU or LSU cross domain),
                    # add the whole group of families like 'LSU,RF02543,RF02540,RF02541,RF02546'
                    for molecule in joint_alignments.keys():
                        if result['accession'] in joint_alignments[molecule]:
                            rfam_families_to_align.add("%s,%s" % (molecule,",".join(joint_alignments[molecule])))

                    #print('Mapping %-12s to %s' % (chain,result['accession']))
                    self.logger.info('Mapping %-12s to %s' % (chain,result['accession']))

            else:
                # record that the chain was considered but not mapped
                if chain in chain_to_range_to_mapping:
                    print('Chain %-12s already mapped to %s, not changing that' % (chain,chain_to_range_to_mapping[chain]['rfam_acc']))
                else:
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

                    chain_to_range_to_mapping[chain] = {}
                    chain_to_range_to_mapping[chain]['        0-        0'] = new_md

                    print('No mapping found for %-12s' % (chain))
                    self.logger.info('No mapping found for %-12s' % (chain))

        # write information about the decisions
        current_datetime = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        with open(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'mapping_decision_%s.txt' % current_datetime),'w') as f:
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
                    rfam_families_to_align.add(chain_to_range_to_mapping[chain][range_string]['rfam_acc']+"_combined")

        # work just with specific Rfam families
        # rfam_families_to_align = set(['RF00023'])  # very slow to align all the sequences!
        rfam_families_to_align = set(['RF01959'])  # good to practice with long alignment
        rfam_families_to_align = set(['RF00005_combined','RF02543_combined'])
        rfam_families_to_align = set(['RF01959_combined'])  # good to practice with long alignment
        rfam_families_to_align = set(['RF02543_combined'])  # uses 50% to 80% memory and then gets killed
        rfam_families_to_align = set(['RF00005_combined'])

        # practice with SSU and LSU alignments
        # rfam_families_to_align = set()
        # for molecule in joint_alignments.keys():
        #     rfam_families_to_align.add("%s,%s" % (molecule,",".join(joint_alignments[molecule])))

        sorted_rfam_families_to_align = sorted(rfam_families_to_align, key = lambda x: ("_" in x,x))

        # write a list of Rfam families that have new chains and so need to be re-aligned
        print('Writing rfam_families_to_align.txt')
        self.logger.info('Writing rfam_families_to_align.txt')
        # only ever append to this list; the only thing that can take away from the list is the aligning program
        with open(rfam_to_align_filename,'w') as f:
            for rfam_family in sorted_rfam_families_to_align:
                if not 'RF00000' in rfam_family and not "_combined" in rfam_family:
                    # write PDB chain sequences to align in this Rfam family
                    print('Writing PDB chain sequences to align to %s' % rfam_family)
                    self.logger.info('Writing PDB chain sequences to align to %s' % rfam_family)
                    num_chains = self.write_pdb_chain_sequences_to_align(RFAM_ALIGNMENT_DIRECTORY,rfam_family,chain_to_range_to_mapping,chain_to_sequence)

                if not "RF00000" in rfam_family:
                    f.write('%s\n' % rfam_family)

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

    def rectify_alignment_file(self,alignment_file,alignment_file_out=None):
        """
        cmalign AFA format splits long sequences onto multiple lines
        Might as well put those on one line here.
        cmalign Pfam format puts header on the same line as sequence; separate those
        Also, replace leading and trailing - with .
        """

        if not alignment_file_out:
            alignment_file_out = alignment_file

        with open(alignment_file,'r') as f:
            lines = f.readlines()

        new_lines = []

        if lines[0].startswith(">"):
            # fasta format
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
            for line in lines:
                # remove trailing \n
                line = line.rstrip('\n')
                if line.startswith("#=GS"):
                    fields = line.split()
                    id = fields[1]
                    header = id + " ".join(fields[3:])
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
            # Pfam format, from how I have heard it described
            k = 0
            i = 0
            while i < len(lines):
                if lines[i].startswith("#") or lines[i].startswith("/") or len(lines[i]) < 5:
                    i += 1
                    continue
                if k == 0:
                    # find the first column of sequence, using the first line of the file
                    k = len(lines[i])-1
                    while k > 0 and lines[i][k] != ' ':
                        k -= 1
                new_lines.append(">"+lines[i][:k].rstrip() + "\n")  # remove spaces after header
                new_sequence = self.rectify_sequence(lines[i][k+1:].replace('\n',''))
                new_lines.append(new_sequence + "\n")
                i += 1

        with open(alignment_file_out,'w') as f:
            f.writelines(new_lines)


    def data(self, rfam_family, **kwargs):
        """
        Input is a list of Rfam families that need to be aligned.
        First plain Rfam families, and we align the PDB chains.
        Then, names like RF00001_combined, and we align PDB chains and sequences.
        That way, the PDB chains get aligned quickly, and then we work on full families, which is slower.
        rfam_family could be a single family or a list for a joint alignment like 'LSU,RF02543,RF02540,RF02541,RF02546'
        """

        # read the current list of Rfam families that need to be aligned
        # if the previous process failed, it will try again now
        rfam_families_to_align = set()
        with open(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'rfam_families_to_align.txt'),'r') as f:
            for line in f.readlines():
                rfam_families_to_align.add(line.replace('\n',''))

        rfam_family_original = rfam_family
        if "_combined" in rfam_family:
            align_combined_file = True
            rfam_family = rfam_family.replace("_combined","")
        else:
            align_combined_file = False

        if 'LSU' in rfam_family:
            model_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'covariance_models','RF02543.cm') # eukaryotic LSU
        elif 'SSU' in rfam_family:
            model_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'covariance_models','RF01960.cm') # eukaryotic SSU
        else:
            model_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'covariance_models','%s.cm' % rfam_family)

        if not os.path.exists(model_file):
            # extract the covariance models that will be needed
            self.extract_covariance_models([rfam_family])

            if not os.path.exists(model_file):
                raise core.Skip("Covariance model for %s is not found" % (rfam_family))

        pdb_chain_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_PDB_chains.fa' % rfam_family.split(',')[0])
        if not os.path.exists(pdb_chain_file):
            raise core.Skip("PDB chain file for %s is not found" % (rfam_family))

        if rfam_family in rfam_families_to_align and not align_combined_file:
            # align PDB chains mapped to this Rfam family
            print("Aligning PDB chains mapped to %s using Infernal cmalign" % rfam_family)

            alignment_file_fa = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_PDB_chains.fa' % rfam_family.split(',')[0])
            alignment_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_PDB_chains.Pfam' % rfam_family.split(',')[0])
            cmalign_location = os.path.join(INFERNAL_LOCATION,'src','cmalign')

            # run infernal; default mxsize is 128, 1200 covers most big alignments
            command = '%s --mxsize 1200 --outformat AFA %s %s > %s' % (cmalign_location,model_file,pdb_chain_file,alignment_file)
            command = '%s --mxsize 1200 --noprob --outformat Pfam %s %s > %s' % (cmalign_location,model_file,pdb_chain_file,alignment_file)
            print('Running command: %s' % command)
            os.system(command)

            # even if Infernal fails, the alignment file will be created, so then a .gz version will be made
            # check to see that Infernal really produced an alignment
            # if so gzip, if not remove
            if os.path.getsize(alignment_file) > 0:
                # read alignment file, put sequence all on one line, replace leading and trailing - with .
                self.rectify_alignment_file(alignment_file,alignment_file_fa)
                os.system("gzip -f %s" % (alignment_file_fa))
                copy_command = "cp %s /var/www/html/data/alignments/rfam/." % (alignment_file + ".gz")
                os.system(copy_command)

                # update the list of rfam families to align to remove this one for PDB chains
                with open(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'rfam_families_to_align.txt'),'w') as f:
                    for rf in sorted(rfam_families_to_align):
                        # remove rfam_family from rfam_families_to_align
                        if not rf == rfam_family:
                            f.write('%s\n' % rf)

            else:
                os.remove(alignment_file)
                self.logger.info('Infernal failed to align %s' % pdb_chain_file)
                self.logger.info('Running command: %s' % command)
                raise core.Skip("Something went wrong with aligning sequences for %s_combined" % (rfam_family))

        elif rfam_family_original in rfam_families_to_align and align_combined_file:

            # align collected sequences to this Rfam family
            # sequence_file_gz = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_rfam_sequences.fa.gz' % rfam_family)

            #bring over previously made sequence files, by removing the PDB chains; no longer needed
            # if not os.path.exists(sequence_file_gz):
            #     old_sequence_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_full_anchored.fa.gz' % rfam_family)
            #     if os.path.exists(old_sequence_file):
            #         with gzip.open(old_sequence_file,'r') as f:
            #             lines = f.readlines()
            #         with gzip.open(sequence_file_gz,'w') as f:
            #             i = 0
            #             found_rfam_sequence = False
            #             while i < len(lines):
            #                 if not found_rfam_sequence:
            #                     fields = lines[i].split(",")[0].split("_")
            #                     if len(fields) == 4:
            #                         # this line and the next are from a PDB chain; skip
            #                         i += 2
            #                     else:
            #                         found_rfam_sequence = True
            #                 if found_rfam_sequence:
            #                     f.write(lines[i])
            #                     f.write(lines[i+1])
            #                     i += 2


            # rewrite the code below to handle the case of multiple Rfam families for SSU and LSU joint alignments

            sequence_file_gz = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_rfam_sequences.fa.gz' % rfam_family)

            if os.path.exists(sequence_file_gz):
                print("Aligning sequences  mapped to %s using Infernal cmalign" % rfam_family)

                # combine PDB chain sequences with Rfam sequences
                combined_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_combined.fa' % rfam_family)

                with open(combined_file,'w') as f:
                    # write pdb chain sequences
                    with open(pdb_chain_file,'r') as pdb_chain_f:
                        for line in pdb_chain_f:
                            f.write(line)
                    # write rfam sequences
                    with gzip.open(sequence_file_gz,'r') as rfam_f:
                        if rfam_family in ['RF02543']:
                            # RF02543 may have a sequence that is too long for Infernal
                            max_length = 0
                            lines = rfam_f.readlines()
                            i = 0
                            while i < len(lines):
                                max_length = max(max_length,len(lines[i+1]))
                                if len(lines[i+1]) < 8000:
                                    f.write(lines[i])
                                    f.write(lines[i+1])
                                i += 2
                            print('Max length of sequence in RF02543 is %i, but only using those up to 8000' % max_length)
                        else:
                            for line in rfam_f:
                                f.write(line)


                # use Pfam format
                alignment_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_combined.Pfam' % rfam_family)
                # use Stockholm format
                alignment_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_combined.sto' % rfam_family)
                # then convert to fasta format in rectify_alignment_file
                alignment_file_fa = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_combined.fa' % rfam_family)
                cmalign_location = os.path.join(INFERNAL_LOCATION,'src','cmalign')

                # run infernal; default mxsize is 128, 1200 covers most big alignments
                command = '%s --mxsize 1200 --noprob %s %s > %s' % (cmalign_location,model_file,combined_file,alignment_file)
                command = '%s --mxsize 1200 --noprob --outformat Pfam %s %s > %s' % (cmalign_location,model_file,combined_file,alignment_file)
                print('Running command: %s' % command)
                os.system(command)

                # even if Infernal fails, the alignment file will be created, so then a .gz version will be made
                # check to see that Infernal really produced an alignment
                # if so gzip, if not remove
                if os.path.getsize(alignment_file) > 0:
                    # read alignment file, put sequence all on one line, replace leading and trailing - with .
                    self.rectify_alignment_file(alignment_file,alignment_file_fa)
                    os.system("gzip -f %s" % (alignment_file_fa))
                    os.system("gzip -f %s" % (combined_file))      # zip the unaligned file of PDB chains and Rfam sequences
                    os.system("gzip -f %s" % (alignment_file))

                    copy_command = "cp %s /var/www/html/data/alignments/rfam/." % (alignment_file_fa + ".gz")
                    os.system(copy_command)

                    # update the list of rfam families to align
                    with open(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'rfam_families_to_align.txt'),'w') as f:
                        for rf in sorted(rfam_families_to_align):
                            # remove rfam_family from rfam_families_to_align
                            if not rf == rfam_family + "_combined":
                                f.write('%s\n' % rf)

                else:
                    os.remove(alignment_file)
                    print('  Infernal failed to align %s' % combined_file)
                    self.logger.info('Running command: %s' % command)
                    raise core.Skip("Something went wrong with aligning sequences for %s_combined" % (rfam_family))

            else:
                # note that the sequence file is needed for this Rfam family
                rfam_sequence_file_needed = set()
                sequence_needed_filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'rfam_sequence_file_needed.txt')
                if os.path.exists(sequence_needed_filename):
                    with open(sequence_needed_filename,'r') as f:
                        for line in f.readlines():
                            rfam_sequence_file_needed.add(line.replace('\n',''))
                rfam_sequence_file_needed.add(rfam_family)
                with open(sequence_needed_filename,'w') as f:
                    for rf in sorted(rfam_sequence_file_needed):
                        f.write('%s\n' % rf)

                print('No sequence file found for %s' % rfam_family)
                raise core.Skip("Rfam sequence file for %s is not found" % (rfam_family))

        else:
            raise core.Skip("Rfam family %s is not in rfam_families_to_align.txt" % (rfam_family))

