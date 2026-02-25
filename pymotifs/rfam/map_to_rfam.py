"""
Get PDB chain sequences that have not been mapped to Rfam.
Run cmsearch to map them.
More.

Rfam.pdb is downloaded from the current version of Rfam
pdb_chain_to_rfam.txt has the same format and includes additional mappings
In both cases, each chain may map to more than one Rfam family
That can happen because a chain has multiple molecules attached into one chain
Also because a chain fits multiple Rfam families reasonably well,
like a bacterial SSU chain might also match the archaeal and eukaryotic SSU families.

A chain that matches multiple Rfam families will be written into multiple alignment files.
So the alignment files are not going to always show the best mapping for a chain,
they will show additional mappings beyond what is stored in the chain_property_value table.
The choice of the best mapping happens in consensus_name.py

"""

from collections import defaultdict
import datetime
import gzip
import os
import requests

import pymotifs.core as core
from pymotifs import models as mod
from sqlalchemy import or_

# constants that should maybe be defined elsewhere
RFAM_ALIGNMENT_DIRECTORY = '/usr/local/pipeline/alignments'
INFERNAL_LOCATION = '/usr/local/pipeline/alignments/infernal-1.1.4'
INFERNAL_LOCATION = '/usr/local/pipeline/alignments/infernal-1.1.5'
MIN_PDB_CHAIN_LENGTH = 20   # this matches what you see in Rfam.pdb

rfam_release = '15.0'

# lists of Rfam families to include in cross-domain alignments
# align according to the covariance model for the first family listed
joint_alignments = {}
joint_alignments['SSU'] = ['RF01960','RF01959','RF00177','RF02545','RF02542'] # align to eukaryotic SSU
joint_alignments['LSU'] = ['RF02543','RF02540','RF02541','RF02546'] # align to eukaryotic LSU

# new in May 2025, let's see how they work
joint_alignments['5.8S-bacteria'] = ['RF02541','RF00002']        # align to bacterial LSU
joint_alignments['5.8S-archaea'] = ['RF02540','RF00002']         # align to archaeal LSU
joint_alignments['LSU-archaea-bacteria'] = ['RF02541','RF02540'] # align to bacterial LSU
joint_alignments['SSU-archaea-bacteria'] = ['RF00177','RF01959'] # align to bacterial SSU

class Loader(core.Loader):
    merge_data = True
    mark = False
    allow_no_data = True

    # need to add stages that produce chain_info
    dependencies = set([])


    def read_mapping(self,filename):
        """
        Read a file like Rfam.pdb or pdb_chain_to_rfam.txt and store the chain to Rfam mappings there
        What is returned is chain_to_range_to_mapping with entries like
            chain_to_range_to_mapping["4V9F|1|0"]["000000001-000002910"] = {}

        Complicating matters, there may be more than one mapping stored per chain
        For example, in 6YDP there are very long chains with multiple matches per chain

        Format:
        rfam_acc	pdb_id	chain	pdb_start	pdb_end	bit_score	evalue_score	cm_start	cm_end	hex_colour
        RF00001	4u3u	7	1	121	88.80	6.5e-20	1	119	8484c0
        """

        if not os.path.exists(filename):
            return {}

        if filename.endswith('.gz'):
            with gzip.open(filename,'rt') as f:
                lines = f.readlines()
        else:
            with open(filename,'rt') as f:
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

                    range_string = '%09s-%09s' % (fields[3],fields[4])  # start and stop positions in the chain

                    if range_string in chain_to_range_to_mapping[chain]:
                        if float(md['bit_score']) > float(chain_to_range_to_mapping[chain][range_string]['bit_score']):
                            # keep the mapping if it has a higher bit score
                            chain_to_range_to_mapping[chain][range_string] = md
                    else:
                        chain_to_range_to_mapping[chain][range_string] = md

                    count += 1

        self.logger.info('Found %d mappings in %s' % (count,filename))

        return chain_to_range_to_mapping

    def write_mapping(self,filename,chain_to_range_to_mapping):
        """
        Write a file with a format like Rfam.pdb

        rfam_acc	pdb_id	chain	pdb_start	pdb_end	bit_score	evalue_score	cm_start	cm_end	hex_colour
        RF00001	4u3u	7	1	121	88.80	6.5e-20	1	119	8484c0
        """

        lines = []
        for chain in chain_to_range_to_mapping.keys():
            for range_string in chain_to_range_to_mapping[chain].keys():
                cm = chain_to_range_to_mapping[chain][range_string]

                header_list = ['rfam_acc','pdb_id','chain','pdb_start','pdb_end','bit_score','evalue_score','cm_start','cm_end','hex_colour']
                line = ''
                for header in header_list:
                    if header in cm:
                        if header == 'pdb_id':
                            line += str(cm[header]).upper()
                        elif header == 'hex_colour' and not cm[header] == 'rna3dhub':
                            # change hex color to rfam to indicate who mapped it
                            line += 'rfam'
                        else:
                            line += str(cm[header])
                    else:
                        print('No data value for header %s in %s' % (header,cm))

                    if header == 'hex_colour':
                        line += '\n'
                    else:
                        line += '\t'
                lines.append(line)

        # sort lines by Rfam family, decreasing score
        lines = sorted(lines, key=lambda x: (x.split("\t")[0],-float(x.split("\t")[5])))
        new_header = "rfam_acc	pdb_id	chain	pdb_start	pdb_end	bit_score	evalue_score	cm_start	cm_end	mapped_by\n"
        lines = [new_header] + lines

        if filename.endswith('.gz'):
            with gzip.open(filename,'w') as f:
                f.writelines(lines)
        else:
            with open(filename,'wt') as f:
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
        family_id = rfam_family.split(",")[0]
        filename = os.path.join(directory,'sequences','%s_PDB_chains.fa' % family_id)

        if num_chains > 0:
            # long headers can make cmalign crash, so write a separate file with just the headers
            header_filename = os.path.join(directory,'sequences','%s_headers.txt' % family_id)
            with open(header_filename,'wt') as f:
                c = 1
                for sequence in sequence_to_chain_range.keys():
                    f.write(">header%d\t>%s\n" % (c,",".join(sorted(sequence_to_chain_range[sequence]))))
                    c += 1

            # write short headers in the sequence file
            with open(filename,'wt') as f:
                c = 1
                for sequence in sequence_to_chain_range.keys():
                    # f.write(">" + ",".join(sorted(sequence_to_chain_range[sequence])) + "\n")
                    f.write(">header%d\n" % c)

                    # replace O with X so that cmalign won't crash
                    sequence = sequence.replace("O","X")
                    f.write(sequence + "\n")
                    c += 1
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
            with gzip.open(filename,'rt') as f:
                lines = f.readlines()
        else:
            with open(filename,'rt') as f:
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
                            if result['score'] >= rfam_family_to_minima[rfam_family]['score'] and score_per_length >= min(0.5,rfam_family_to_minima[rfam_family]['score_per_length']):
                                # new result is better than existing result
                                results.append(result)
                        else:
                            # new assignments not present in Rfam, be conservative
                            if result['score'] > 30 and score_per_length > 0.3 and length > MIN_PDB_CHAIN_LENGTH:
                                results.append(result)
                            # code below resulted in poor matches with low score and short match length
                            # elif score_per_length > 0.8:
                            #     results.append(result)

        # go from best match (highest bit score) to worst match, keep non-overlapping matches for each chain
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


    def map_chains_to_best_mapping(self, chain_to_range_to_mapping):
        """
        Choose the best mapping for each chain, when there are multiple choices.
        Return chain_to_range_to_best_mapping where chain is like "4V9F|1|0",
        range is like 000000004-000000056 and the value is a dictionary
        each entry of chain_to_range_to_mapping[chain][range] is a dictionary with fields
        rfam_acc	pdb_id	chain	pdb_start	pdb_end	bit_score	evalue_score	cm_start	cm_end	hex_colour
        """

        chain_to_range_to_best_mapping = {}
        for chain in chain_to_range_to_mapping:
            # Keep the best non-overlapping ranges for this chain
            # get all range,mapping tuples
            range_mapping = chain_to_range_to_mapping[chain].items()
            # sort ranges from best bit score to worst
            range_mapping = sorted(range_mapping, key=lambda x: float(x[1]['bit_score']), reverse=True)
            indices_seen = set()
            keep_range_mapping = []
            for r,mapping in range_mapping:
                if not mapping["rfam_acc"] == "RF00000":
                    pdb_start = int(mapping["pdb_start"])
                    pdb_end   = int(mapping["pdb_end"])
                    new_indices = set(range(pdb_start,pdb_end+1))

                    # allow small overlaps between ranges
                    if len(new_indices & indices_seen) < 10:
                        # have not seen this range before
                        keep_range_mapping.append((r,mapping))
                        indices_seen = indices_seen | new_indices

            if len(keep_range_mapping) > 0:
                chain_to_range_to_best_mapping[chain] = {}
                if len(keep_range_mapping) == 1:
                    r, mapping = keep_range_mapping[0]
                    chain_to_range_to_best_mapping[chain][r] = mapping
                else:
                    self.logger.info('Multiple ranges in %s' % chain)
                    # query database to find chain indices with resolved nucleotides
                    fields = chain.split("|")
                    pdb_id = fields[0]
                    chain_name = fields[2]
                    with self.session() as session:
                        query = session.query(mod.UnitInfo.chain_index).\
                            filter(mod.UnitInfo.pdb_id == pdb_id).\
                            filter(mod.UnitInfo.chain == chain_name)
                        observed_indices = set([row.chain_index for row in query])

                    for r,mapping in keep_range_mapping:
                        pdb_start = int(mapping["pdb_start"])
                        pdb_end   = int(mapping["pdb_end"])
                        observed_indices_mapped = observed_indices & set(range(pdb_start,(pdb_end+1)))
                        self.logger.info("%4d observed nts in mapping %s" % (len(observed_indices_mapped),mapping))

                        if len(observed_indices_mapped) > 10:
                            # enough observed nucleotides to keep this chain to range to mapping
                            chain_to_range_to_best_mapping[chain][r] = mapping

        return chain_to_range_to_best_mapping


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
            with open(rfam_to_align_filename,'rt') as f:
                for line in f.readlines():
                    rfam_families_to_align.add(line.replace('\n',''))

        # read Rfam.pdb to identify chains that Rfam has mapped
        # Rfam.pdb needs to be downloaded with each new Rfam release
        self.logger.info('Reading mappings from Rfam.pdb')
        chain_to_range_to_mapping = self.read_mapping(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'Rfam.pdb'))

        # find minimum score in each Rfam family according to Rfam mapping
        self.logger.info('Getting minimum Rfam scores in each family')
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

        self.logger.info('Rfam.pdb has %i chains from %d families' % (len(chain_to_range_to_mapping.keys()),len(rfam_family_to_minima.keys())))

        # read pdb_chain_to_rfam.txt to identify chains that are already mapped by Rfam or rna3dhub
        # multiple mappings of the same chain may be present
        # if the same range is mapped multiple times, the one with the higher bit score is used
        self.logger.info('Reading previous mappings')
        rna3dhub_chain_to_range_to_mapping = self.read_mapping(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'pdb_chain_to_rfam.txt'))

        self.logger.info('pdb_chain_to_rfam.txt has %i chains' % len(rna3dhub_chain_to_range_to_mapping.keys()))

        # if a chain is not mapped by Rfam but was mapped here, add it
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
                                  mod.ChainInfo.compound,
                                  mod.PdbInfo.title,
                                  mod.PdbInfo.release_date).\
                outerjoin(mod.PdbInfo, mod.ChainInfo.pdb_id == mod.PdbInfo.pdb_id).\
                filter(or_(mod.ChainInfo.entity_macromolecule_type.contains('olyribonucleotide'),
                           mod.ChainInfo.entity_macromolecule_type.contains('ybrid')))

        chain_to_sequence = {}
        sequence_to_chains = defaultdict(list)  # many chains have the same sequence
        chain_to_title = {}
        chain_to_compound = {}
        chain_to_release_date = {}

        for r in query:
            chain = "%s|1|%s" % (r.pdb_id.upper(),r.chain_name)
            chain_to_sequence[chain] = r.sequence
            chain_to_title[chain] = r.title
            chain_to_compound[chain] = r.compound
            chain_to_release_date[chain] = str(r.release_date)

            if len(r.sequence) > MIN_PDB_CHAIN_LENGTH:
                sequence_to_chains[r.sequence].append(chain)

        # if any unmapped sequences are the same as a chain that is already mapped, map them like before
        # this avoids having to run cmsearch on them, and promotes consistency.  Also propagates any errors!
        self.logger.info('Identifying chains whose sequences are already mapped to Rfam')
        for sequence, chains in sequence_to_chains.items():
            new_chains = set(chains) - set(chain_to_range_to_mapping.keys())
            mapped_chains = set(chains) & set(chain_to_range_to_mapping.keys())

            # don't bother checking, because we know the same SSU sequence can map to multiple
            # SSU families with different ranges, nothing more controversial is likely to happen
            # if len(mapped_chains) > 1:
            #     # it's conceivable that the same sequence has different mapping information
            #     # do some simple checks to see if they are the same
            #     range_to_count = defaultdict(int)
            #     for chain in mapped_chains:
            #         for range_string in chain_to_range_to_mapping[chain].keys():
            #             range_to_count[range_string] += 1
            #     counts = list(range_to_count.values())
            #     if min(counts) < max(counts):
            #         self.logger.info('Chains %s have the same sequence but different mappings' % sorted(mapped_chains))
            #         print('Chains %s have the same sequence but different mappings' % sorted(mapped_chains))

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
                            cm[key] = value
                        cm['pdb_id'] = chain.split('|')[0]
                        cm['chain'] = chain.split('|')[2]
                        chain_to_range_to_mapping[chain][range_string] = cm

                    rfam_mapped_to = [chain_to_range_to_mapping[chain][range_string]['rfam_acc'] for range_string in chain_to_range_to_mapping[chain].keys()]
                    print('Chain %s has same sequence as %s which maps to %s' % (chain,mapped_chain,rfam_mapped_to))
                    self.logger.info('Chain %s has same sequence as %s which maps to %s' % (chain,mapped_chain,rfam_mapped_to))

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
            with open(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'chains_for_cmsearch.fa'),'wt') as f:
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
        with open(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'mapping_decision','mapping_decision_%s.txt' % current_datetime),'wt') as f:
            # write basic information about each Rfam family's mappings
            for rfam_family, min_score in sorted(rfam_family_to_minima.items()):
                f.write('%s Min score %7.2f Min length %5d Min ratio %7.2f\n' % (rfam_family,min_score['score'],min_score['length'],min_score['score_per_length']))
            for output_line in sorted(output_lines):
                f.write(output_line)

        # write all data to pdb_chain_to_rfam.txt
        # may include multiple mappings for some chains, especially those directly from Rfam
        self.logger.info("Writing any new lines to pdb_chain_to_rfam.txt")
        print('Writing mappings to pdb_chain_to_rfam.txt')
        self.logger.info('chain_to_range_to_mapping has %i chains' % len(chain_to_range_to_mapping.keys()))
        self.write_mapping(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'pdb_chain_to_rfam.txt'),chain_to_range_to_mapping)

        # read previous file of best mappings; these are best and also approved mappings
        previous_chain_to_range_to_best_mapping = self.read_mapping(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'pdb_chain_to_best_rfam.txt'))
        previous_rfam_to_chains = {}
        for chain in previous_chain_to_range_to_best_mapping.keys():
            for r in previous_chain_to_range_to_best_mapping[chain].keys():
                rfam = previous_chain_to_range_to_best_mapping[chain][r]['rfam_acc']
                if not rfam in previous_rfam_to_chains:
                    previous_rfam_to_chains[rfam] = set()
                previous_rfam_to_chains[rfam].add(chain)

        # identify the best non-overlapping ranges for each chain
        # also verify that there are observed nucleotides in that range, to avoid cases where
        # the experimental sequence contains molecules that do not have x,y,z coordinates
        chain_to_range_to_best_mapping = self.map_chains_to_best_mapping(chain_to_range_to_mapping)

        # identify current mapping of rfam to chains
        rfam_to_chains = {}
        for chain in chain_to_range_to_best_mapping.keys():
            for r in chain_to_range_to_best_mapping[chain].keys():
                rfam = chain_to_range_to_best_mapping[chain][r]['rfam_acc']
                if not rfam in rfam_to_chains:
                    rfam_to_chains[rfam] = set()
                rfam_to_chains[rfam].add(chain)

        # identify Rfam families that have at least one new chain, and so need to be re-aligned
        for rfam in rfam_to_chains.keys():
            if not rfam in previous_rfam_to_chains or not previous_rfam_to_chains[rfam] == rfam_to_chains[rfam]:
                rfam_families_to_align.add(rfam)

        # extend Rfam families with new chain to joint alignments containing that Rfam family
        # add the whole group of families like 'LSU,RF02543,RF02540,RF02541,RF02546'
        for molecule in joint_alignments.keys():
            if set(joint_alignments[molecule]) & rfam_families_to_align:
                rfam_families_to_align.add("%s,%s" % (molecule,",".join(joint_alignments[molecule])))

        # when you need to re-align everything, you can use this list
        # most of the time this should be "if False:"
        if True:
            rfam_families_to_align = set()
            for chain in chain_to_range_to_best_mapping.keys():
                for range_string in chain_to_range_to_best_mapping[chain].keys():
                    rfam_families_to_align.add(chain_to_range_to_best_mapping[chain][range_string]['rfam_acc'])

        # work just with specific Rfam families
        # rfam_families_to_align = set(['RF00386'])
        # rfam_families_to_align = set(['RF00005','RF02543'])
        # rfam_families_to_align = set(['RF03064','RF03022','RF03072','RF00061','RF03125','RF03121'])

        # re-make joint alignments
        # for molecule,families in joint_alignments.items():
        #     rfam_families_to_align.add("%s,%s" % (molecule,",".join(families)))

        rfam_families_to_align = sorted(rfam_families_to_align)

        # write a list of Rfam families that have new chains and so need to be re-aligned
        print('Writing rfam_families_to_align.txt')
        self.logger.info('Writing rfam_families_to_align.txt')
        # only ever append to this list here
        # the only thing that can take away from the list is the aligning program
        with open(rfam_to_align_filename,'wt') as f:
            for rfam_family in list(rfam_families_to_align):
                if not 'RF00000' in rfam_family:
                    # write PDB chain sequences to align in this Rfam family
                    num_chains = self.write_pdb_chain_sequences_to_align(RFAM_ALIGNMENT_DIRECTORY,rfam_family,chain_to_range_to_best_mapping,chain_to_sequence)
                    if num_chains > 0:
                        print('Wrote %4d PDB chain sequences to align to %s' % (num_chains,rfam_family))
                        self.logger.info('Wrote %4d PDB chain sequences to align to %s' % (num_chains,rfam_family))
                        f.write('%s\n' % rfam_family)

        # identify Rfam families that are manually approved to have mappings
        with open('/usr/local/pipeline/hub-core/aux/manual_consensus_names.tsv','rt') as f:
            lines = f.readlines()
        approved_rfam = set()
        for line in lines:
            rfam = line.strip().split("\t")[0]
            if rfam.startswith('RF'):
                approved_rfam.add(rfam)

        # families that are mapped by rna3dhub are tracked in the file rfam_mapped_by_rna3dhub.txt
        # with open('/usr/local/pipeline/hub-core/aux/rfam_mapped_by_rna3dhub.txt','rt') as f:
        #     lines = f.readlines()
        # mapped_by_rna3dhub = set()
        # for line in lines:
        #     rfam = line.strip()
        #     if rfam.startswith('RF'):
        #         mapped_by_rna3dhub.add(rfam)

        # identify families that are mapped by RNA3DHub but not by Rfam
        mapped_by_rfam = set(rfam_family_to_minima.keys())
        mapped_by_rna3dhub = approved_rfam - mapped_by_rfam
        with open('/usr/local/pipeline/hub-core/aux/rfam_mapped_by_rna3dhub.txt','wt') as f:
            f.writelines("\n".join(sorted(mapped_by_rna3dhub)))

        # make note of any families that are not properly accounted for
        # mapped_by_rna3dhub_but_not_noted = approved_rfam - mapped_by_rna3dhub - mapped_by_rfam
        # self.logger.info('mapped_by_rna3dhub_but_not_noted %s' % sorted(mapped_by_rna3dhub_but_not_noted))

        # mapped_by_rna3dhub_but_not_approved = mapped_by_rna3dhub - approved_rfam
        # self.logger.info('mapped_by_rna3dhub_but_not_approved %s' % sorted(mapped_by_rna3dhub_but_not_approved))

        mapped_by_rfam_but_no_standardized_name = mapped_by_rfam - approved_rfam
        self.logger.info('mapped_by_rfam_but_no_standardized_name %s' % mapped_by_rfam_but_no_standardized_name)
        with open('/usr/local/pipeline/hub-core/aux/mapped_by_rfam_but_no_standardized_name.txt','wt') as f:
            f.writelines("\n".join(sorted(mapped_by_rfam_but_no_standardized_name)))

        with open('/usr/local/pipeline/hub-core/aux/rfam_mapped_by_rfam.txt','wt') as f:
            f.writelines("\n".join(sorted(mapped_by_rfam)))

        # identify potential mappings that should be checked manually
        # cull out the mappings that are not manually approved
        mappings_to_check_manually = []
        for chain in list(chain_to_range_to_best_mapping.keys()):
            for r in list(chain_to_range_to_best_mapping[chain].keys()):
                mapping = chain_to_range_to_best_mapping[chain][r]
                if not mapping['rfam_acc'] in approved_rfam:
                    mappings_to_check_manually.append(mapping)
                    # remove this mapping from the output so we have a solid list of
                    # best and approved mappings for JAR3D and other purposes
                    del chain_to_range_to_best_mapping[chain][r]

        if len(mappings_to_check_manually) > 0:
            # get a mapping of Rfam accession ids to Rfam family name
            # RF00001 -> 5S ribosomal RNA
            rfam_family_file = "/usr/local/pipeline/alignments/family.txt"
            with open(rfam_family_file,'rt') as f:
                lines = f.readlines()
            rfam_id_to_name = {}
            for line in lines:
                fields = line.split('\t')
                name = fields[3]
                rfam_family = fields[0]
                rfam_id_to_name[rfam_family] = name

            # write a file of potential mappings that should be checked manually
            mappings_to_check_manually = sorted(mappings_to_check_manually, key = lambda x : (x['rfam_acc'],x['pdb_id']))
            header = "pdb_id\tchain\tPDB start\tPDB end\tbit_score\tcm start\tcm end\tRfam ID\tRfam name\tPDB chain compound name\tStandardized name(s)\tDomain; subgroup\tPDB release date\tPDB structure title"
            with open("/usr/local/pipeline/hub-core/aux/rfam_check_manually.tsv","wt") as f:
                f.writelines(header+"\n")
                for new_data in mappings_to_check_manually:
                    # rfam_acc	pdb_id	chain	pdb_start	pdb_end	bit_score	evalue_score	cm_start	cm_end	hex_colour
                    chain = "%s|1|%s" % (new_data['pdb_id'].upper(),new_data['chain'])
                    line  = new_data['pdb_id'] + '\t'
                    line += new_data['chain'] + '\t'
                    line += new_data['pdb_start'] + '\t'
                    line += new_data['pdb_end'] + '\t'
                    line += "%0.2f\t" % float(new_data['bit_score'])
                    line += new_data['cm_start'] + '\t'
                    line += new_data['cm_end'] + '\t'
                    line += new_data['rfam_acc'] + '\t'
                    line += rfam_id_to_name[new_data['rfam_acc']] + '\t'
                    line += chain_to_compound[chain] + '\t\t\t'        # leave empty columns
                    line += chain_to_release_date[chain] + '\t'
                    line += chain_to_title[chain]

                    f.writelines(line+"\n")

        # write list of best mappings; these also satisfy the minimum requirements for mappings
        # these are also approved mappings
        self.logger.info("Writing best mapping from chain to Rfam family")
        print("Writing best mapping from chain to Rfam family")
        self.write_mapping(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'pdb_chain_to_best_rfam.txt'),chain_to_range_to_best_mapping)

        if len(rfam_families_to_align) == 0:
            print('No new Rfam families to align')
            raise core.Skip('No new Rfam families to align')

        return rfam_families_to_align

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

        with gzip.open(rfam_model_file,'rt') as file:
            cm_data = file.read().split("//")

        for i in range(len(cm_data)):
            model = cm_data[i].strip()
            if model.startswith("INFERNAL"):
                for family in rfam_families:
                    if family in model:
                        family_model_filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'covariance_models','%s.cm' % family)
                        with open(family_model_filename, 'wt') as out_file:
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


    def read_stockholm_file(self,alignment_file):
        id_to_header = {}
        id_to_sequence = {}
        id_list_sequence = []
        id_list_pdb = []

        with open(alignment_file,'rt') as f:
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

        return id_to_header, id_to_sequence, id_list_sequence, id_list_pdb


    def read_fasta_file(self,fasta_file):
        with gzip.open(fasta_file,'rt') as f:
            lines = f.readlines()

        id_to_fasta_sequence = {}

        for this_line in lines:
            line = this_line.strip("\n")
            if line.startswith(">"):
                id = line[1:]
            else:
                id_to_fasta_sequence[id] = line

        return id_to_fasta_sequence


    def map_sequence_position_to_column(self,line,starting_position=1,ending_position=float('Inf')):
        """
        Internally, column counting starts at 0 because these are indexes into Python strings
        01234
        .......AUGCUC-GGGUAAUCGCUGCGGCCGGu--------------------------------------uUCGG-----CCGUA

        starting_position tells what sequence position in the PDB sequence is the first
        one in the alignment
        """

        gap_characters = '-.'
        column_to_sequence_position = {}
        sequence_position_to_column = {}

        column_print_limit = 0
        if column_print_limit > 0:
            print(line[:column_print_limit])

        sequence_position = starting_position
        for column in range(len(line)):
            if not line[column] in gap_characters and sequence_position <= ending_position:
                # first sequence_position is 1 because that is how RNA3DHub lists them
                column_to_sequence_position[column] = sequence_position
                sequence_position_to_column[sequence_position] = column

                if column < column_print_limit:
                    print("%4d %4d %s" % (column,sequence_position,line[column]))

                sequence_position += 1

        return column_to_sequence_position, sequence_position_to_column


    def map_rfam_seed_columns_to_unit_ids(self,rfam_family,source='rfam'):
        """
        Needs to work either with original Rfam seed alignments or with CaCoFold seed alignments.
        Logic is the same after that.

        """

        if source == 'rfam':
            rfam_seed_alignment_filename = os.path.join("/usr/local/pipeline/alignments/alignments",rfam_family+"_rfam_seed.sto")
        else:
            rfam_seed_alignment_filename = os.path.join("/usr/local/pipeline/alignments/alignments",rfam_family+"_cacofold.sto")

        if not os.path.exists(rfam_seed_alignment_filename):
            self.logger.info('Cannot find %s' % rfam_seed_alignment_filename)
            return

        if source == 'rfam':
            combined_alignment_filename  = os.path.join("/usr/local/pipeline/alignments/alignments",rfam_family+"_combined_seed.fa.gz")
        else:
            combined_alignment_filename  = os.path.join("/usr/local/pipeline/alignments/alignments",rfam_family+"_combined_cacofold.fa.gz")

        if not os.path.exists(combined_alignment_filename):
            self.logger.info('Cannot find %s' % combined_alignment_filename)
            return

        id_to_header, id_to_sequence, id_list_sequence, id_list_pdb = self.read_stockholm_file(rfam_seed_alignment_filename)

        id_to_fasta_sequence = self.read_fasta_file(combined_alignment_filename)

        # map seed column to combined alignment column, sequence by sequence
        seed_column_to_combined_column = {}
        for id in id_to_sequence.keys():
            # print(id_to_sequence[id])
            # print(id_to_fasta_sequence.get(id,"no sequence"))

            # work with the seed sequence line
            # loop over seed column number (starting at 0, internally)
            # map column (starting at 0) to sequence position (starting at 1)
            # map sequence position (starting at 1) to alignment column (starting at 0)
            seed_line = id_to_sequence[id]
            seed_column_to_sequence_position, seed_sequence_position_to_column = self.map_sequence_position_to_column(seed_line)

            # work with the fasta sequence line
            # loop over alignment column number (starting at 0, internally)
            # map sequence position (starting at 1) to alignment column (starting at 0)
            combined_line = id_to_fasta_sequence[id]
            combined_column_to_sequence_position, combined_sequence_position_to_column = self.map_sequence_position_to_column(combined_line)

            for seed_column, seed_sequence_position in seed_column_to_sequence_position.items():
                combined_column = combined_sequence_position_to_column[seed_sequence_position]
                if not seed_column in seed_column_to_combined_column:
                    seed_column_to_combined_column[seed_column] = combined_column
                else:
                    if not combined_column == seed_column_to_combined_column[seed_column]:
                        print('Disagreement about seed column %d' % seed_column)
                        input("Press Enter")

        # for seed_column, combined_column in sorted(seed_column_to_combined_column.items()):
        #     print('Seed column %3d corresponds to combined column %3d in %s' % (seed_column,combined_column,rfam_family))

        # read pdb_chain_to_best_rfam.txt to know what PDB chain sequence position the alignment starts at
        chain_to_range_to_best_mapping = self.read_mapping(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'pdb_chain_to_best_rfam.txt'))

        # loop over 3D structure chains, map seed column to combined column to sequence position to unit id
        seed_column_to_unit_ids = {}
        chain_list = []
        for chains, sequence in id_to_fasta_sequence.items():
            if not chains in id_to_sequence:
                # this line of the combined alingnment is for one or more PDB chains
                for chain_id in chains.split(","):
                    # download sequence position to unit id mappings for this chain
                    fields = chain_id.split("_")
                    pdb = fields[0]
                    chain = fields[1]
                    chain_triple = "%s|1|%s" % (pdb.upper(),chain)

                    self.logger.info('rfam_family %s' % rfam_family)
                    self.logger.info('chain_triple %s' % chain_triple)

                    # get starting sequence position of this chain in the alignment
                    starting_position = None
                    ending_position = None
                    for range, mapping_dictionary in chain_to_range_to_best_mapping.get(chain_triple,{}).items():
                        if mapping_dictionary['rfam_acc'] == rfam_family:
                            fields = range.split("-")
                            if len(fields) == 2:
                                starting_position = int(fields[0])
                                ending_position = int(fields[1])

                    if starting_position and ending_position:
                        chain_list.append(chain_triple)

                        # self.logger.info('range and dictionary %s' % chain_to_range_to_best_mapping.get(chain_triple,{}))
                        # self.logger.info('start %s end %s' % (starting_position,ending_position))

                        # process the PDB sequence on this line of the combined alignment
                        combined_column_to_sequence_position, combined_sequence_position_to_column = self.map_sequence_position_to_column(sequence,starting_position,ending_position)

                        url = "https://rna.bgsu.edu/rna3dhub/rest/SeqtoUnitMapping?ife=%s" % (chain_triple)
                        response = requests.get(url)
                        sequence_position_to_unit_id = {}
                        if response.status_code == 200:
                            # data is tab-delimited lines
                            for line in response.text.split("</br>"):
                                # print(line)
                                fields = line.split()
                                if len(fields) == 3:
                                    sequence_fields = fields[0].split("|")
                                    sequence_position = int(sequence_fields[4])  # starting at 1
                                    unit_id = fields[2]
                                    # if unit_id == "NULL":
                                    #     unit_id = ""
                                    if sequence_position in sequence_position_to_unit_id:
                                        if len(unit_id) < len(sequence_position_to_unit_id[sequence_position]):
                                            # use the shorter unit id
                                            sequence_position_to_unit_id[sequence_position] = unit_id
                                    else:
                                        sequence_position_to_unit_id[sequence_position] = unit_id

                        # loop over seed columns, map as far as possible
                        # every seed column needs an entry for every 3D structure to keep unit ids aligned
                        for seed_column, combined_column in seed_column_to_combined_column.items():
                            unit_id = ''
                            if combined_column in combined_column_to_sequence_position:
                                sequence_position = combined_column_to_sequence_position[combined_column]
                                if sequence_position in sequence_position_to_unit_id:
                                    unit_id = sequence_position_to_unit_id[sequence_position]
                                    # print('Seed column %d to combined column %d to sequence position %d to unit id %s in %s' % (seed_column,combined_column,sequence_position,unit_id,rfam_family))
                            if not seed_column in seed_column_to_unit_ids:
                                seed_column_to_unit_ids[seed_column] = []
                            seed_column_to_unit_ids[seed_column].append(unit_id)

        # write out a file mapping seed column to unit ids
        if source == 'rfam':
            output_filename  = os.path.join("/usr/local/pipeline/alignments/alignments","%s_%s_column_to_unit_id.txt" % (rfam_family,rfam_release))
        else:
            output_filename  = os.path.join("/usr/local/pipeline/alignments/alignments","%s_%s_CCF_column_to_unit_id.txt" % (rfam_family,rfam_release))

        with open(output_filename,'wt') as f:
            t = "%s\t%s\n" % (rfam_family,"\t".join(chain_list))
            f.writelines(t)
            for seed_column in sorted(seed_column_to_unit_ids.keys()):
                unit_ids = seed_column_to_unit_ids[seed_column]
                # output file seed alignment column numbering starts at 1 to match JAR3D and normal counting
                t = "%d\t%s\n" % (seed_column+1,"\t".join(unit_ids))
                f.writelines(t)


    def rectify_alignment_file(self,alignment_file,header_filename,alignment_file_out_gz=None):
        """
        Read alignment in various formats, write out in fasta format
        cmalign sto format splits long sequences onto multiple lines; combine
        cmalign AFA format splits long sequences onto multiple lines; combine
        cmalign Pfam format has long header on one line, then short header and sequence later
        Replace leading and trailing - with .
        Replace short placeholder headers with full header of PDB chain list
        """

        if not alignment_file_out_gz:
            alignment_file_out_gz = alignment_file+".gz"

        short_header_to_full = {}
        with open(header_filename,'rt') as f:
            for line in f:
                fields = line.split()
                if len(fields) == 2:
                    short_header = fields[0]
                    full_header = fields[1]
                    short_header_to_full[short_header] = full_header

        new_lines = []
        if ".fa" in alignment_file.lower():
            # fasta format, not so bad to read all of the lines at once
            with open(alignment_file,'rt') as f:
                lines = f.readlines()

            i = 0
            while i < len(lines):
                if lines[i].startswith(">"):
                    if lines[i] in short_header_to_full:
                        new_lines.append(short_header_to_full[lines[i]])
                    else:
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

            id_to_header, id_to_sequence, id_list_sequence, id_list_pdb = self.read_stockholm_file(alignment_file)

            for id in id_list_pdb + id_list_sequence:
                if ">"+id_to_header[id] in short_header_to_full:
                    new_lines.append(short_header_to_full[">"+id_to_header[id]] + "\n")
                else:
                    new_lines.append(">%s\n" % id_to_header[id])
                sequence = id_to_sequence[id]
                new_sequence = self.rectify_sequence(sequence)
                new_lines.append(new_sequence + "\n")

        elif ".pfam" in alignment_file.lower():
            # Pfam format, from how I have heard it described, pretty much Stockholm
            k = 0
            with open(alignment_file,'rt') as f:
                for line in f:
                    if line.startswith("#") or line.startswith("/") or len(line) < 5:
                        continue
                    if k == 0:
                        # find the first column of sequence, using the first line of the file
                        k = len(line)-1
                        while k > 0 and line[k] != ' ':
                            k -= 1
                    header = ">"+line[:k].rstrip()    # remove spaces after header
                    if header in short_header_to_full:
                        new_lines.append(short_header_to_full[header] + "\n")
                    else:
                        new_lines.append(header + "\n")
                    new_sequence = self.rectify_sequence(line[k+1:].replace('\n',''))
                    new_lines.append(new_sequence + "\n")

        with gzip.open(alignment_file_out_gz,'wt') as f:
            f.writelines(new_lines)


    def align_progressive(self,rfam_family,model_file,sequence_file,alignment_file_sto):
        """
        Align sequences in this family progressively, a set number at a time.
        Also useful for aligning some files of PDB chains.
        sequence_file is the sequences to align
        alignment_file_sto is the name of the output file
        """

        cmalign_location = os.path.join(INFERNAL_LOCATION,'src','cmalign')
        alimerge_location = os.path.join(INFERNAL_LOCATION,'easel','miniapps','esl-alimerge')

        # read the entire sequence file
        if sequence_file.endswith('.gz'):
            with gzip.open(sequence_file,'rt') as f:
                lines = f.readlines()
        else:
            with open(sequence_file,'rt') as f:
                lines = f.readlines()

        self.logger.info('Found %i lines in large sequence family %s' % (len(lines),sequence_file))
        print('Found %i lines in large sequence family %s' % (len(lines),sequence_file))

        # lines_delta must be even!  because the code counts lines, not header-sequence units
        if rfam_family == 'RF00005':
            lines_delta = 100000
        elif rfam_family == 'RF02543':
            lines_delta = 200
        elif rfam_family == 'LSU':
            lines_delta = 50         # 25 sequences at a time
        elif rfam_family == 'SSU':
            lines_delta = 50         # 25 sequences at a time
        else:
            lines_delta = 2*int(len(lines)/20.0)

        # small temporary file for alignment of a portion of the sequences
        sequence_temp = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_progressive_temp.fa' % rfam_family)

        # growing alignment file
        alignment_growing = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_progressive_growing.sto' % rfam_family)

        # merged alignment file
        alignment_merged = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_progressive_merged.sto' % rfam_family)

        start_index = 0
        while start_index < len(lines):
            last_index = min(start_index+lines_delta,len(lines))
            if len(lines) - last_index < lines_delta*0.4:
                # just do all the rest of the sequences
                last_index = len(lines)
                lines_delta = last_index - start_index

            with open(sequence_temp,'wt') as f:
                f.writelines(lines[start_index:last_index])

            # first time write to the final output location "merged", after that write to the "temp" file
            if start_index == 0:
                alignment_temp = alignment_merged
            else:
                alignment_temp = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_progressive_temp.sto' % rfam_family)

            # run infernal; default mxsize is 128, 1200 covers most big alignments
            # mxsize 1000 worked better for LSU_PDB_chains on 2024-07-08
            # create alignment_temp file
            if rfam_family in ['RF02541']:
                # not sure why a medium-length sequence requires 1260 mb but let's try this
                alignment_command = '%s --mxsize 1300 --noprob --outformat Pfam %s %s > %s' % (cmalign_location,model_file,sequence_temp,alignment_temp)
            else:
                alignment_command = '%s --mxsize 1000 --noprob --outformat Pfam %s %s > %s' % (cmalign_location,model_file,sequence_temp,alignment_temp)
            self.logger.info('Start %5d last %5d running command: %s' % (start_index,last_index,alignment_command))
            print('Start line %5d last %5d running command: %s' % (start_index,last_index,alignment_command))
            os.system(alignment_command)

            # check to see that Infernal really produced an alignment
            # if not, remove and bail out
            if os.path.exists(alignment_temp) and os.path.getsize(alignment_temp) == 0:
                os.remove(alignment_temp)
                self.logger.info('  Infernal did not align %s' % sequence_temp)
                print('  Infernal did not align %s' % sequence_temp)
                self.logger.info('  Failed command was: %s' % command)
                raise core.Skip("Something went wrong with aligning sequences for %s; try again next time" % (rfam_family))

            # merge the temporary new alignment with the growing full alignment
            if start_index > 0:
                # move the previously produced large alignment so it can be added to
                command = 'mv %s %s' % (alignment_merged,alignment_growing)
                os.system(command)

                # merge alignment_temp and alignment_growing to produce new alignment_merged
                command = '%s --outformat pfam -o %s %s %s' % (alimerge_location,alignment_merged,alignment_temp,alignment_growing)
                print('Combining alignments by running command: %s' % command)
                os.system(command)

            start_index += lines_delta

        # if this completed successfully, move the merged file to the final location
        if os.path.exists(alignment_merged):
            if os.path.exists(alignment_file_sto):
                os.remove(alignment_file_sto)
            command = 'mv %s %s' % (alignment_merged,alignment_file_sto)
            os.system(command)

        os.remove(alignment_temp)
        os.remove(alignment_growing)

        return alignment_command


    def data(self, rfam_family, **kwargs):
        """
        Input is a single Rfam family that needs to be aligned.
        First we align the PDB chains.
        That way, the PDB chains get aligned quickly, and then we work on full families, which is slower.
        rfam_family could be a single family or a list for a joint alignment like 'LSU,RF02543,RF02540,RF02541,RF02546'
        """

        rfam_fields = rfam_family.split(",")

        family_id = rfam_fields[0]  # e.g., LSU, RF02543

        pdb_alignment_file_sto = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_PDB_chains.sto' % family_id)
        alimerge_location = os.path.join(INFERNAL_LOCATION,'easel','miniapps','esl-alimerge')

        align_pdb_chains = True

        if align_pdb_chains:
            if len(rfam_fields) > 1:
                # first listed family is the covariance model to align to
                model_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'covariance_models','%s.cm' % (rfam_fields[1]))
            else:
                model_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'covariance_models','%s.cm' % rfam_family)

            cmalign_location = os.path.join(INFERNAL_LOCATION,'src','cmalign')

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
            # LSU pdb chains killed on 2025-05-21
            # SSU pdb chains segmentation fault on 2025-05-21

            if family_id in ['LSU','SSU']:
                # use progressive alignment
                command = self.align_progressive(family_id,model_file,pdb_chain_file,pdb_alignment_file_sto)
            else:
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
            self.logger.info('Infernal not able to align %s' % pdb_chain_file)
            self.logger.info('Command was: %s' % command)
            raise core.Skip("Something went wrong with aligning PDB chain sequences for %s" % (rfam_family))
        else:
            # read PDB chain alignment file, put sequence all on one line, replace leading and trailing - with .
            # output as .fa.gz file
            alignment_file_fa = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_PDB_chains.fa' % family_id)
            header_filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_headers.txt' % family_id)
            self.rectify_alignment_file(pdb_alignment_file_sto,header_filename,alignment_file_fa+".gz")

            combined_sto = ""
            if len(rfam_fields) == 1:
                # combine PDB chain alignment with seed alignments from Rfam, CaCoFold
                # for source in ['rfam','cacofold']:
                for source in ['cacofold']:
                    # check to see if Rfam seed alignment file is downloaded
                    if source == 'rfam':
                        seed_alignment_file_sto = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_rfam_seed.sto' % rfam_family)
                        if not os.path.exists(seed_alignment_file_sto):
                            self.logger.info("Downloading Rfam seed alignment for %s" % rfam_family)
                            print("Downloading Rfam seed alignment for %s" % rfam_family)

                            url = "https://rfam.org/family/%s/alignment/stockholm" % rfam_family

                            response = requests.get(url)

                            if response.status_code == 200:
                                with open(seed_alignment_file_sto, "wt") as file:
                                    file.write(response.text)
                            else:
                                self.logger.info('Could not download from %s' % url)

                    else:
                        seed_alignment_file_sto = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_cacofold.sto' % rfam_family)
                        # if we don't have the file, there is no way to recover here

                    if os.path.exists(seed_alignment_file_sto):
                        # merge PDB chain alignment with Rfam seed alignment
                        if source == 'rfam':
                            combined_seed_sto = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_combined_seed.sto' % family_id)
                        else:
                            combined_seed_sto = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_combined_cacofold.sto' % family_id)

                        command = '%s --outformat pfam -o %s %s %s' % (alimerge_location,combined_seed_sto,pdb_alignment_file_sto,seed_alignment_file_sto)
                        self.logger.info('Combining PDB chain and Rfam seed alignments')
                        print('Combining PDB chain and Rfam seed alignments')
                        print('Running command: %s' % command)
                        os.system(command)

                        # read combined alignment file, put sequence all on one line, replace leading and trailing - with .
                        self.logger.info("Rectifying alignment of PDB chains and seed")
                        print("Rectifying alignment of PDB chains and seed")
                        if source == 'rfam':
                            seed_alignment_file_fa = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_combined_seed.fa' % family_id)
                        else:
                            seed_alignment_file_fa = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_combined_cacofold.fa' % family_id)
                        header_filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_headers.txt' % family_id)
                        self.rectify_alignment_file(combined_seed_sto,header_filename,seed_alignment_file_fa+".gz")

                        # map columns of the Rfam seed alignment to unit ids
                        self.map_rfam_seed_columns_to_unit_ids(rfam_family,source)

                if False:

                    # check to see if sequences from Rfam full family are already collected; if not, note that
                    sequence_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_rfam_sequences.fa' % rfam_family)
                    sequence_file_gz = sequence_file + ".gz"

                    if not os.path.exists(sequence_file) and not os.path.exists(sequence_file_gz):
                        # note that the sequence file is needed for this Rfam family
                        rfam_sequence_file_needed = set()
                        sequence_needed_filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'rfam_sequence_file_needed.txt')
                        # read list of Rfam families that need sequence files
                        if os.path.exists(sequence_needed_filename):
                            with open(sequence_needed_filename,'rt') as f:
                                for line in f.readlines():
                                    rfam_sequence_file_needed.add(line.replace('\n',''))
                        rfam_sequence_file_needed.add(rfam_family)
                        # write list of Rfam families that need sequence files
                        with open(sequence_needed_filename,'wt') as f:
                            for rf in sorted(rfam_sequence_file_needed):
                                f.write('%s\n' % rf)

                        # self.logger.info('No sequence file %s found for %s' % (sequence_file,rfam_family))
                        print('No sequence file found for %s' % rfam_family)
                        raise core.Skip("Rfam sequence file %s for %s is not found; check again on the next run" % (sequence_file,rfam_family))

                    # check to see if sequences from Rfam full family are already aligned; if not, align them
                    seq_alignment_file_sto = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_rfam_sequences.sto' % rfam_family)
                    seq_alignment_file_sto_gz = seq_alignment_file_sto + ".gz"

                    if not os.path.exists(seq_alignment_file_sto) and not os.path.exists(seq_alignment_file_sto_gz):
                        self.logger.info("Aligning sequences from %s full family using Infernal cmalign" % rfam_family)
                        print("Aligning sequences from %s full family using Infernal cmalign" % rfam_family)

                        # some Rfam families have so many sequences, or such long sequences, that they need a different approach
                        if rfam_family in ['RF00005','RF02541','RF02543']:
                            model_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'covariance_models','%s.cm' % rfam_family)
                            sequence_file = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_rfam_sequences.fa.gz' % rfam_family)
                            self.align_progressive(rfam_family,model_file,sequence_file,seq_alignment_file_sto)

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
                    self.logger.info("Rectifying full alignment")
                    print("Rectifying alignment")
                    alignment_file_fa = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'alignments','%s_combined.fa' % family_id)
                    header_filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,'sequences','%s_headers.txt' % family_id)
                    self.rectify_alignment_file(combined_sto,header_filename,alignment_file_fa+".gz")

            # read the current list of Rfam families that need to be aligned
            rfam_families_to_align = set()
            with open(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'rfam_families_to_align.txt'),'rt') as f:
                for line in f.readlines():
                    rfam_families_to_align.add(line.replace('\n',''))

            # remove rfam_family from the list of rfam families to align, if everything above was successful
            with open(os.path.join(RFAM_ALIGNMENT_DIRECTORY,'rfam_families_to_align.txt'),'wt') as f:
                for rf in sorted(rfam_families_to_align):
                    # remove rfam_family from rfam_families_to_align
                    if not rf == rfam_family:
                        f.write('%s\n' % rf)

            # clean up files we don't need
            if len(combined_sto) > 0:
                os.remove(combined_sto)
            os.remove(pdb_alignment_file_sto)

            if False and len(rfam_fields) == 1:
                os.system("gzip -f %s" % (seq_alignment_file_sto))

