"""The common exceptions that show up across the pipeline.
"""


class StageFailed(Exception):
    """This is raised when one stage of the pipeline fails.
    """
    def __init__(self, msg, ids=None):
        super(StageFailed, self).__init__(msg)
        self.failed_ids = ids


class InvalidState(Exception):
    """This is an exception meant to be used when we have entered into some
    sort of invalid state in a stage. For example if some data which is
    required to be in the database is not present then this should be raised.
    Or when the stage should produce data, but nothing is saved then this is
    raised.
    """
    pass


class Skip(Exception):
    """Class to indicate that the processing of one input in for one stage
    should be skipped.
    """
    pass
