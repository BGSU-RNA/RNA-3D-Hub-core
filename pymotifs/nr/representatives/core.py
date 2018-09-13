import abc
import copy
import operator as op

from pymotifs import core

from pymotifs.constants import NR_ALLOWED_METHODS


class Representative(core.Base):
    """Just a base class for all things that find representatives. This is only
    for book keeping purposes and implements nothing.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractproperty
    def method(self):
        """
        The method name for this representative finder. This should be a
        unique string among all representative subclasses.
        """
        return None

    def filter_group_by_method(self, group, methods=NR_ALLOWED_METHODS):
        """
        This will filter the group that is being examiend to just a copy of
        the parent, and members entries. This is done so that we have a group
        that we can manipulate without messing up parts elsewhere. In addition,
        the members of the group will be filtered to only those with allowed
        methods. These are the methods listed in the methods set.

        Parameters
        ----------
        group : dict
            A group dictonary that must have 'parent' and 'members' entries.
        methods : set
            A set of method names that are allowed. If no members have the
            given method then all are used.

        Returns
        -------
        copied : dict
            A group dictonary with only the 'parent' and 'members' entries.
        """

        meth = op.itemgetter('method')
        members = [ife for ife in group['members'] if meth(ife) in methods]
        if not members:
            members = group['members']

        return {
            'parent': copy.deepcopy(group.get('parent', [])),
            'members': copy.deepcopy(members),
        }

    def insert_as_representative(self, representative, members, sort=None):
        """
        This will insert the selected representative at the head of the list
        of members. If sort is given it will also sort the remaining members
        using the given function.

        Parameters
        ----------
        representative : object
            The selected representative
        members : list
            List of other members of the class
        sort : lambda, None
            The function to use to sort the members if needed.
        """

        ordered = [representative]
        to_add = [m for m in members if m['id'] != representative['id']]
        if sort:
            to_add.sort(key=sort, reverse=True)
        ordered.extend(to_add)
        return ordered
