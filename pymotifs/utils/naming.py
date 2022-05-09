import random
import itertools as it

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
            self.logger.debug("Group with %i members", len(group['members']))

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
            self.logger.info("Named group %s with %i members" %
                (named_group['comment'], len(named_group['members'])))

            named.append(named_group)
            handles.add(named_group['name']['handle'])

        if len(named) != len(groups):
            raise core.InvalidState("Missing groups in naming")

        return named


class ChangeCounter(core.Base):
    """A class to help summarize changes between two sets of named groups.
    """

    def group_changes(self, groups, parent_groups):
        """Compute the number of changes at the group level. That is the number
        of added, removed, updated or unchanged groups.
        """

        as_handle = lambda g: g['name']['handle']
        as_name = lambda g: (g['name']['handle'], g['name']['version'])
        parents = set(as_name(g) for g in parent_groups)
        parent_handles = set(as_handle(g) for g in parent_groups)

        added = []
        updated = []
        unchanged = []
        for group in groups:
            name = as_name(group)
            if name in parents:
                unchanged.append(group)
            elif group['name']['handle'] in parent_handles:
                updated.append(group)
            else:
                added.append(group)

        removed = []
        handles = set(as_handle(g) for g in groups)
        for parent in parent_groups:
            handle = as_handle(parent)
            if handle not in handles:
                removed.append(parent)

        return {
            'added': added,
            'removed': removed,
            'updated': updated,
            'unchanged': unchanged
        }

    def transformed_changes(self, groups, parents, fn):
        """Compute the counts of changes given some transformation function.

        :param list groups: The list of groups to compare.
        :param list parents: The list of parent groups to compare against.
        :param function fn: The function to use for transforming.
        :returns: A dictonary of the
        """

        def as_set(entries):
            mapped = it.imap(fn, entries)
            return set(it.chain.from_iterable(mapped))

        current = as_set(groups)
        parent = as_set(parents)

        return {
            'added': current - parent,
            'removed': parent - current,
            'unchanged': current.intersection(parent),
        }

    def __as_counts__(self, changes):
        counts = {}
        for name, change in changes.items():
            counts[name] = len(change)
        return counts

    def __call__(self, groups, parents, **transformers):
        """Compute the number of changes between groups and parent groups.
        """

        def members(group):
            return [m['id'] for m in group['members']]

        group_changes = self.group_changes(groups, parents)
        member_changes = self.transformed_changes(groups, parents, members)

        data = {
            'groups': self.__as_counts__(group_changes),
            'members': self.__as_counts__(member_changes),
        }

        keys = ['added', 'unchanged', 'updated']
        total = sum(len(group_changes[n]) for n in keys)
        assert len(groups) == total

        for entry_name, fn in transformers.items():
            changes = self.transformed_changes(groups, parents, fn)
            data[entry_name] = self.__as_counts__(changes)

        return data
