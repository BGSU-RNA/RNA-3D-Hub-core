import csv

from pymotifs import core
from pymotifs import models as mod
from sqlalchemy import desc

class Loader(core.SimpleLoader):

    merge_data = True
    allow_no_data = True

    """
    @property
    def table(self):
        return mod.LoopAnnotations
    """

    def to_process(self, pdbs, **kwrags):
        """ return the different loop types to process, see if that works
        """
        return ["IL"]

    """
    def remove(self, *args, **kwargs):
        Does not actually remove from the DB. We always want the loop ids to
        be consistent so we do not automatically remove loops.
        """

    """Added has_data because it should force the stage to recompute "IL". This is not a permenant solution."""
    
    def has_data(self, args, **kwargs):
        return False
    
    def query(self, session, pdb):
        return session.query(mod.LoopAnnotations)

    def _get_loop_annotations(self):
        """ Gets access to the spreadsheet, reads values, and saves data in DB friendly format.
        """
        
        """Create dictionary for loop_annotations 1 and 2 and motif_annotations 1 and 2 with loop_id and motif_id 
        as their respective keys and with values consisting of a 2-element set: annotation and author.
        """

        loop_annotation1 = {}
        loop_annotation2 = {}
        motif_annotation1 = {}
        motif_annotation2 = {}

        #Read data from csv file and save each row as a list of strings

        print("annotations.py reading loop annotations")

        with open("/usr/local/pipeline/hub-core/IL_annotations.csv", "r") as annotations:
            entries = csv.reader(annotations, delimiter = ',')
            for row in entries:
                #If a loop_id is indicated, then check for and save loop annotation 1 and/or loop annotation 2
                if isinstance(row[1], str) and row[1]:
                    if isinstance(row[2], str) and row[2]:
                        loop_annotation1[row[1]]= [row[2], row[4]]
                    if isinstance(row[3], str) and row[3]:
                        loop_annotation2[row[1]]= [row[3], row[4]]
                        
                #Otherwise save motif annotation 1 and/or motif annotation 2       
                else:
                    if isinstance(row[2], str) and row[2]:
                        motif_annotation1[row[0]] = [row[2], row[4]]
                        
                    if isinstance(row[3], str) and row[3]:
                        motif_ananotation2[row[0]] = [row[3], row[4]]

        #Save each loop annnotation to loop_data in DB-friendly format
        
        loop_data = []

        #if a loop instance has both annotation1 and annotation 2 indicated, then do not check for motif annotations
        both_annotations = set(loop_annotation1.keys().intersection(loop_annotation2.keys()))
        
        for loop in both_annotations:
            loop_data.append(mod.LoopAnnotations(
                loop_id = loop,
                annotation_1 = loop_annotation1[loop][0],
                annotation_2 = loop_annotation2[loop][0],
                author = loop_annotation1[loop][1]))

        #The rest of the LoopAnnotations entries need to be checked for motif_annotations
        
        #Get current motif release
        with self.session() as session:
            current_ml_release = session.query(mod.MlReleases.ml_release_id).\
                order_by(desc(mod.MlReleases.date)).\
                limit(1)

        print(current_ml_release)


        #For each motif, get the associated loop_ids.
        all_motifs = set(motif_annotation1.keys().union(motif_annotation2.keys()))
        for motif in all_motifs:
            with self.session() as session:
                query = session.query(mod.MlLoops.loop_id).\
                     filter(mod.MlLoops.motif_id == motif).\
                     filter(mod.MlLoops.ml_release_id == current_ml_release)
                loop_ids = [r.loop_id for r in query]
            
            #For each loop in this motif - check for loop instance annotations, otherwise use the motif annotation (if it exists)
            #This code should handle the cases that the loop has either annotation1, annotation2, or neither
            for loop in loop_ids:
                #If a loop in this motif does not have a loop instance annotation1, then use motif_annotation1 given that motif_annotation1 exists
                if loop not in loop_annotation1.keys():
                    if motif in motif_annotation1:
                        loop_annotation1[loop] = motif_annotation1[motif]
                    else:
                        loop_annotation1[loop] = ["NULL", "NULL"]
                #If a loop in this motif does not have a loop instance annotation2, then use motif_annotation2 given that motif_annotation2 exists
                if loop not in loop_annotation2.keys(): 
                    if motif in motif_annotation2:
                        loop_annotation2[loop] = motif_annotation2[motif]
                    else:
                        loop_annotation2[loop] = ["NULL", "NULL"]

                #Add completed loop_entry to loop_data
                loop_data.append(mod.LoopAnnotations(
                    loop_id = loop,
                    annotation_1 = loop_annotation1[loop][0],
                    annotation_2 = loop_annotation2[loop][0],
                    author = loop_annotation1[loop][1]))
             
        return loop_data
    
    def data(self, pdb, **kwargs):
        """ Get the loop annotations from the spreadsheet.
        Return a list of all loop annotations. """

        print("annotations.py running data method")

        #Get loop annotations
        data = []
        data = self._get_loop_annotations()
#Return data to be written to the database.
        if not data:
            self.logger.error("No loop annotations found.")
        else:
            return data
