import random

from pymotifs import core

from pymotifs.models import NrClasses


class Namer(core.Base):
    size_cutoff = 2.0 / 3.0

    def overlap(self, group1, group2):
        """Compute the overlap between the two groups. This is done by
        comparing the number of autonomous unit ids in common between them
        """

        autonomous1 = set(mem['id'] for mem in group1['members'])
        autonomous2 = set(mem['id'] for mem in group2['members'])

        intersection = autonomous1.intersection(autonomous2)
        if intersection:
            return {'group': group2, 'intersection': intersection}

        return {}

    def overlap_size(self, group, overlap):
        inter_size = len(overlap['intersection'])
        size = max(len(group['members']), len(overlap['group']['members']))
        return float(inter_size) / float(size)

    def same_name(self, group):
        return {
            'handle': group['handle'],
            'version': group['version'],
            'comment': 'Exact match'
        }

    def updated_name(self, group, parents):
        parent_comment = '1 parent'
        if parents == 2:
            parent_comment = '2 parents'

        return {
            'handle': group['handle'],
            'version': group['version'] + 1,
            'comment': 'Updated, %s' % parent_comment
        }

    def new_name(self, count):
        """Create a new name for a class. This will generate a new handle that
        has never been seen before and create a dictonary with the parts named.

        :returns: A naming dictonary.
        """

        known = set()
        with self.session() as session:
            query = session.query(NrClasses.handle).distinct()
            known = set(res.handle for res in query)

        handle = None
        while not handle or handle not in known:
            handle = '%05d' % random.randrange(99999)

        if not count:
            parent_comment = 'no parents'
        elif count == 1:
            parent_comment = '1 parent'
        elif count == 2:
            parent_comment = '%i parents' % count
        else:
            parent_comment = '> 2 parents'

        return {
            'handle': handle,
            'version': 1,
            'comment': 'New id, %s' % parent_comment
        }

    def one_parent(self, group, parent):
        # - use same name if the groups are identical
        if self.overlap_size(group, parent) == 1:
            print(parent)
            return self.same_name(parent['group']['name'])

        # - use a new name if the overlap is > 2/3
        elif self.change_size(group, parent) >= self.size_cutoff:
            return self.updated_name(parent['group']['name'], 1)

        return self.new_name(1)

    def two_parents(self, group, parents):
        parent = min(parents, key=self.overlap_size)
        if self.overlap_size(group, parent) >= self.size_cutoff:
            return self.updated_name(parent, 2)
        return self.new_name(2)

    def parents(self, group, known_groups):
        parents = []
        for known in known_groups:
            overlap = self.overlap(known, group)
            if overlap:
                parents.append(overlap)
        return parents

    def __call__(self, groups, known_groups):
        named = []
        for group in groups:
            parents = self.parents(group, known_groups)

            # No overlaps means new group thus new name
            name = {}
            if not parents:
                name = self.new_name(0)

            elif len(parents) == 1:
                name = self.one_parent(group, parents[0])

            elif len(parents) == 2:
                name = self.two_parents(group, parents)

            # If there is more than 2 parents we always use a new name
            else:
                name = self.new_name(len(parents))

            named_group = dict(group)
            named_group['parents'] = parents
            named_group.update(name)
            named.append(named_group)

        return named
