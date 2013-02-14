"""

Due to a bug, loop HL_4HXH_001 didn't match itself in a pairwise search
and an entry was added to a corresponding no_candidates.txt file.

This script scans and corrects no_candidates.txt files.

"""

import os
import sys
import glob
import pdb
import shutil


# path to all-against-all searches
aAa = ''

# pdb ids of files for which aAa searches should be purged
bad = []

for loop_dir in glob.glob(aAa + '*'):
    loop_id = os.path.basename(loop_dir)

    for pdb_id in bad:
        if pdb_id in loop_id:
            print loop_id
            no_candidates = os.path.join(aAa, loop_id, 'no_candidates.txt')
            # copy for backup
#             shutil.copy(no_candidates, '/Users/anton/Desktop/backup_bad_no_candidates/' + loop_id + '_no_candidates.txt')
            # erase
            f = open(no_candidates, 'w')
            f.close()
            break

sys.exit()

to_fix = []

for loop_dir in glob.glob(aAa + '*'):
    loop_id = os.path.basename(loop_dir)
    no_candidates = os.path.join(aAa, loop_id, 'no_candidates.txt')

    if os.path.exists(no_candidates):
        ids = open(no_candidates).read().splitlines()
        # check if all entries are unique
        uniques = set(ids)
        if len(uniques) != len(ids):
            print '%s: %s vs %s entries' % (loop_id, len(uniques), len(ids))
        # check if its own loop id is in the no_candidates list
        if loop_id in ids:
            print 'Incorrect entry found %s' % loop_id
            to_fix.append(loop_id)

            f = open(no_candidates,'w')
            ids.remove(loop_id)
            for loop in ids:
                f.write(loop + "\n")
            f.close()

    else:
        print 'No_candidates.txt not found for %s' % loop_id

print to_fix