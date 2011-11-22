"""Motif Loader
Loads RNA motif clusters into a database. Manages id tracking, version numbers etc.

Usage: python aPyMotifLoader [options]

Options:
  -f ..., --file=...    use specified csv file with motifs
  -n                    increment the release counter. If omitted - the same
                        release number will be used. 
                        1.0, 1.1 vs. 2.0, 2.1 etc
  -m ...                specify short version description
  -h, --help            show this help  

Examples:
python aPyMotifLoader.py -f motifs.csv -n -m 'first release'
"""

import pdb
import sys
import getopt

import aMLDbConnector as adb # anton database
import aMLRelease as amr # anton motif release

__author__ = 'Anton Petrov'



def load_new_release(file, newRelease, message):
    """
    """
    
    db = adb.aDbConnector()	
    
    oldRelease = db.get_last_release()

    newRelease = amr.Release()
    newRelease.load_from_file(file)    
    newRelease.set_version(oldRelease.version + 1)
    newRelease.set_release(newRelease)
    newRelease.set_message(message)
    
    pdb.set_trace()
    if oldRelease.version == 0:

        db.upload_first_release(newRelease)    
    
    else:
        
        mergedRelease = compare_releases(oldRelease, newRelease, db)
        db.load_release(mergedRelease)

                
    db.close() 
    print 'Done'



def main(argv):
    """
    """
    try:
        opts, args = getopt.getopt(argv, "f:m:n", ['help','file='])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
            
    newRelease = False
    message = ''        
    file = ''   
    
    for opt, arg in opts:
        if opt == '-f':
            file = arg
        elif opt in ('-h', '--help'):
            usage()
            sys.exit()
        elif opt == '-m':
            message = arg
        elif opt == '-n':
            newRelease = True

    if file == '':
        usage()
        sys.exit(2)
    else:
        load_new_release(file, newRelease, message)


def usage():
    print __doc__


if __name__ == "__main__":
    main(sys.argv[1:])


    