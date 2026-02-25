"""
Read manual annotations into the loop_annotations table in the database.
Annotations are specified by loop or by motif group.
Precedence is given to loop-specific annotations and to more recent annotations.
"""

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

    def to_process(self, pdbs, **kwargs):
        """
        return the different loop types to process
        """
        return ["HL","IL","JL"]

    """
    def remove(self, *args, **kwargs):
        Does not actually remove from the DB. We always want the loop ids to
        be consistent so we do not automatically remove loops.
        """

    """Added has_data because it should force the stage to recompute "IL". This is not a permanant solution."""
    def has_data(self, args, **kwargs):
        return False

    def query(self, session, pdb):
        return session.query(mod.LoopAnnotations)

    def data(self, loop_type, **kwargs):
        """
        Read CSV file for loop_type.
        Create dictionary for loop_annotations 1 and 2 and motif_annotations 1 and 2 with loop_id and motif_id
        as their respective keys and with values consisting of a 2-element set: annotation and author.
        Return a list of all loop annotations in format the database will save.
        """

        loop_to_annotation1 = {}
        loop_to_annotation2 = {}
        motif_annotation1 = {}
        motif_annotation2 = {}

        # Read data from csv file and save each row as a list of strings

        self.logger.info("Reading %s annotations from CSV file" % loop_type)
        print("annotations.py reading %s loop annotations" % loop_type)

        with open("/usr/local/pipeline/hub-core/%s_annotations.tsv" % loop_type, "r") as annotations:
            # entries = csv.reader(annotations, delimiter = ',')
            raw_entries = annotations.readlines()

        entries = []
        for line in raw_entries:
            entry = line.rstrip("\n").split("\t")
            entries.append(entry)

        # ['IL_48459.1', 'IL_4R0D_018', 'Isolated non-canonical cWW pair', '', 'sree', '2021-03-23']
        # sort entries so ones without specified loops come in the first group,
        # and entries with specified loops come in the second group
        # within the group that has specific loops, sort by increasing date in column 5
        # so later dates overwrite earlier ones
        # within the group that does not have specific loops, sort by decreasing date,
        # so that later dates do not overwrite earlier ones

        specific_loops = []
        motif_groups = []
        for e in entries:
            if len(e[1]) > 5:
                specific_loops.append(e)
            else:
                motif_groups.append(e)
        specific_loops.sort(key=lambda x: x[5])
        motif_groups.sort(key=lambda x: x[5], reverse=True)
        sorted_entries = motif_groups + specific_loops

        # process all entries, with later entries overwriting earlier ones
        changed_group_annotations = []
        for entry in sorted_entries:
            self.logger.info("Reading entry %s" % entry)

            # trim leading and trailing white space
            row = [r.strip() for r in entry]

            # record loop-specific annotations if present
            if isinstance(row[1], str) and row[1]:
                loop_ids = row[1].replace(" ","").replace(";",",").split(",")
                for loop_id in loop_ids:
                    if isinstance(row[2], str) and row[2]:
                        loop_to_annotation1[loop_id] = [row[2], row[4]]     # save annotation and author
                        if isinstance(row[3], str) and row[3]:
                            loop_to_annotation2[loop_id] = [row[3], row[4]] # save annotation and author
                        else:
                            loop_to_annotation2[loop_id] = [row[2], row[4]] # repeat annotation1

                    # Listing a loop id with two blank annotations means
                    # not to annotate that loop, even with the motif annotation
                    if len(row[2]) == 0 and len(row[3]) == 0:
                        loop_to_annotation1[loop_id] = ["NULL", row[4]]
                        loop_to_annotation2[loop_id] = ["NULL", row[4]]

            # Otherwise check for motif group and save motif annotation 1 and/or motif annotation 2
            elif isinstance(row[0],str) and row[0]:
                motif_id = row[0].replace(" ","").replace(";",",")
                if isinstance(row[2], str) and row[2]:
                    motif_annotation1[motif_id] = [row[2], row[4]]   # save annotation and author
                if isinstance(row[3], str) and row[3]:
                    motif_annotation2[motif_id] = [row[3], row[4]]  # save annotation and author

                with self.session() as session:
                    query = session.query(mod.MlLoops.loop_id).\
                        filter(mod.MlLoops.motif_id == motif_id)
                    loop_ids = list(set([r.loop_id for r in query]))

                self.logger.info("Unique loop ids in motif group %s are:" % motif_id)
                self.logger.info(loop_ids)

                # For each loop_id in this motif, save the annotations
                # This code should handle the cases that the loop has either annotation1, annotation2, or neither
                for loop_id in loop_ids:
                    if motif_id in motif_annotation1:
                        if loop_id in loop_to_annotation1:
                            if not loop_to_annotation1[loop_id][0] == motif_annotation1[motif_id][0]:
                                # previous group annotation is being overwritten by a new group annotation
                                # that may not be desirable
                                changed_group_annotations.append('Changed group annotation for loop %s in %s from %s to %s' % (loop_id, motif_id, loop_to_annotation1[loop_id][0], motif_annotation1[motif_id][0]))
                        loop_to_annotation1[loop_id] = motif_annotation1[motif_id]
                        self.logger.info("Loop %s has annotation 1 %s" % (loop_id, motif_annotation1[motif_id]))
                    else:
                        loop_to_annotation1[loop_id] = ["NULL", "NULL"]
                        self.logger.info("Loop %s has annotation 1 NULL" % (loop_id))
                    if motif_id in motif_annotation2:
                        loop_to_annotation2[loop_id] = motif_annotation2[motif_id]
                    else:
                        loop_to_annotation2[loop_id] = ["NULL", "NULL"]

        # if a loop instance has either annotation1 or annotation 2, write to the database
        annotated_loops = list(set(loop_to_annotation1.keys()) | set(loop_to_annotation2.keys()))
        loop_data = []
        for loop_id in annotated_loops:
            if loop_id in loop_to_annotation1:
                ann1 = loop_to_annotation1[loop_id][0]
                author = loop_to_annotation1[loop_id][1]
            else:
                ann1 = "NULL"

            if loop_id in loop_to_annotation2:
                ann2 = loop_to_annotation2[loop_id][0]
                author = loop_to_annotation1[loop_id][1]
            else:
                ann2 = "NULL"

            if not (ann1 == "NULL" and ann2 == "NULL"):
                # at least one non-trivial annotation
                loop_data.append(mod.LoopAnnotations(
                    loop_id = loop_id,
                    annotation_1 = ann1,
                    annotation_2 = ann2,
                    author = author))

        for c in changed_group_annotations:
            self.logger.info(c)

        self.logger.info('Returning %d loop annotations' % len(loop_data))
        print("annotations.py returning %d loop annotations" % len(loop_data))

        return loop_data
