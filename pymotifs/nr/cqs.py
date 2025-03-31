"""
This stage saves composite quality score and other data to the
table nr_cqs.
The data is computed earlier in using_quality and cached.
"""

from pymotifs import core
from pymotifs import models as mod
from pymotifs.constants import NR_CACHE_NAME

from pymotifs.nr.utils import BaseLoader
from pymotifs.nr.classes import Loader as NrClassLoader
from pymotifs.nr.parent_counts import Loader as CountLoader


class Loader(BaseLoader):
    """
    Loader to store quality data for an input equivalence class
    in table nr_cqs.
    """

    dependencies = set([NrClassLoader,CountLoader])

    """
    We allow this to merge data since sometimes we want to replace.
    """
    merge_data = True

    """We allow for no data to be written when appropriate"""
    allow_no_data = True
    mark = False

    def has_data(self, *args, **kwargs):
        """
        We always want stages.py to think that data needs to be computed
        for this stage, so it does not get skipped.
        This works better than trying to manipulate the query method.
        """
        return False


    def data(self, release, **kwargs):
        """
        Collect composite quality scoring data for the equivalence class.

        Parameters
        ----------
        nr_name : str
            Name of class for which to collect NR-level composite
            quality score (CQS) data.

        Returns
        -------
            The required data for the database update step.
        """

        # retrieve the current data on equivalence classes
        data = self.cached(NR_CACHE_NAME)

        grouping = data['groups']

        # loop over equivalence classes
        for group in grouping:
            if group.get('new_nr_cqs_entry',False):
                # loop over ifes in the equivalence class
                for member in group['members']:
                    yield mod.NrCqs(
                        ife_id = member['id'],
                        nr_name = group['name']['full'],
                        maximum_experimental_length = member['max_observed'],
                        fraction_unobserved = member['fraction_unobserved'],
                        percent_observed = member['percent_observed'],
                        composite_quality_score = member['composite_quality_score'])

        # the field name maximum_experimental_length is poorly chosen,
        # because the value is actually the maximum number of observed nucleotides
        # in the unit_info table, not the length of the longest experimental
        # sequence that went into the experiment
        # Maximum number of observed nucleotides is more robust, because sometimes
        # a very long chain gets reported as the experimental sequence, and that
        # includes many chains concatenated together.
