# Possible changes

Here is a list of some possible improvements. They are in no particular order.

1. Set all NULL f_crossing to 0 or recompute them
    NULL in f_crossing should be 0. It may be best to recompute the data with
    NULL's though.

2. Add foreign key from unit_pairs_interactions.f_lwbp to bp_family_info
    bp_family_info is where all bp's are defined so it makes sense to have an
    FK to that table.

3. Prevent NULL in f_crossing
    It should never be null as it is defined for all interactions.

4. Add a unique constraint in exp_seq_unit_mapping (unit_id, exp_seq_position_id)
    Each unit id may only be part of 1 experimental sequence position so it is
    useful to enforce this in the database.

5. Improve normalization procedures. 1OB5 is listed as having no valid
   sequences but it does. This is because the sequence ends FC.

6. Parallelize running of PDB's within a stage. Most stages (Loader,
   SimpleLoader) are parallelizble as each individual PDB can be run without
   dependencies on any other PDB. This should give us a large speed up.

7. Parallelize the running of some stages at the same level. The toposort of
   all stages produces stages in a series of levels. Within each level there
   are no dependencies among the stages so they can be run independently.
   Running them in parallel would speed things up.

8. Split unit_pairs_interactions into several tables. We could have one table
   for base pairs, one for stacks, etc. This would allow us to add a NOT NULL
   constraint to interaction columns. Also, the tables wouldn't be filled with
   columns that have lots of NULL's as they are currently. A view could be made
   to preserve the previous querying method.

9. Cleanup and simplify the hierarchy in pymotifs.core.stages. Currently, there
   are too many levels of inheritance. This can be cut down quite a bit using
   composition instead.

10. Do a better job at not recomputing stages which produce no data. For
    example, units.incomplete will often produce nothing for a given structure
    (eg. 3BO1) and thus each time the pipeline is run the stage will attempt to
    recompute missing for the file. This slows things down as we are redoing
    work that will not produce anything. This can probably be over come by
    using the marks in pdb_analysis_status to determine if this has been
    processed before.

11. Alter the logic that determines what needs to be run so that it does 2
    things, first it checks what is is done and then runs the normal
    to_process. Then it gets the difference between the two and runs those.
    This would remove the need for a 'has_data' method as we know that each
    entry being processed will never have data.

12. Unify all configuration data. There are 3 places where things are
    configured, the configuration file (`conf/motifatlas.json`), the
    pymotifs.constants file, options given on the command line, the values set
    using '--manual key=value' on the command line. It would be good if all
    three were a single thing. Right now each totally separate but it would be
    really useful if any could be set in the same way.
