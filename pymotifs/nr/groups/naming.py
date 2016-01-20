import random

from pymotifs import core


class Namer(core.Base):
    size_cutoff = 2.0 / 3.0

    def overlap(self, group1, group2):
        """Compute the overlap between the two groups. This is done by
        comparing the number of ife ids in common between them.
        """

        members1 = set(mem['id'] for mem in group1['members'])
        members2 = set(mem['id'] for mem in group2['members'])

        intersection = members1.intersection(members2)
        if not intersection:
            return {}

        return {'group': group2, 'intersection': intersection}

    def overlap_size(self, group, overlap):
        inter_size = len(overlap['intersection'])
        size = max(len(group['members']), len(overlap['group']['members']))
        return float(inter_size) / float(size)

    def same_name(self, group):
        return {
            'handle': group['handle'],
            'version': group['version'],
            'comment': 'Exact match',
            'type': 'exact',
        }

    def updated_name(self, group, parents):
        parent_comment = '1 parent'
        if parents == 2:
            parent_comment = '2 parents'

        return {
            'handle': group['handle'],
            'version': group['version'] + 1,
            'comment': 'Updated, %s' % parent_comment,
            'type': 'updated',
        }

    def new_name(self, count, known):
        """Create a new name for a class. This will generate a new handle that
        has never been seen before and create a dictonary with the parts named.

        :returns: A naming dictonary.
        """

        handle = '%05d' % random.randrange(99999)
        while handle in known:
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
            'comment': 'New id, %s' % parent_comment,
            'type': 'new',
        }

    def one_parent(self, group, parent, known):
        # - use same name if the groups are identical
        if self.overlap_size(group, parent) == 1:
            return self.same_name(parent['group']['name'])

        # - use a new name if the overlap is > 2/3
        elif self.overlap_size(group, parent) >= self.size_cutoff:
            return self.updated_name(parent['group']['name'], 1)

        return self.new_name(1, known)

    def two_parents(self, group, parents, known):
        parent = max(parents, key=lambda p: self.overlap_size(group, p))
        if self.overlap_size(group, parent) >= self.size_cutoff:
            return self.updated_name(parent['group']['name'], 2)
        return self.new_name(2, known)

    def many_parents(self, group, parents, known):
        # If there is more than 2 parents we always use a new name
        return self.new_name(len(parents), known)

    def parents(self, group, known_groups):
        parents = []
        for known in known_groups:
            overlap = self.overlap(group, known)
            if overlap:
                parents.append(overlap)
        return parents

    def __call__(self, groups, parent_groups, handles):
        named = []
        for group in groups:
            parents = self.parents(group, parent_groups)
            self.logger.info("Group with %i members", len(group['members']))

            # No overlaps means new group thus new name
            name = {}
            if not parents:
                name = self.new_name(0, handles)

            elif len(parents) == 1:
                name = self.one_parent(group, parents[0], handles)

            elif len(parents) == 2:
                name = self.two_parents(group, parents, handles)

            else:
                name = self.many_parents(group, parents, handles)

            named_group = dict(group)
            named_group['parents'] = [p['group'] for p in parents]
            named_group['comment'] = name.pop('comment')
            named_group['name'] = dict(name)
            self.logger.info("Named group with %i members", len(named_group['members']))

            named.append(named_group)
            handles.add(named_group['name']['handle'])

        if len(named) != len(groups):
            raise core.InvalidState("Missing groups in naming")

        return named
