"""A module to determine which stages of the pipeline to run. The Dispatcher
class is what will determine each stage to run and then run them in the correct
order.
"""

import logging
import itertools as it

from pymotifs import core
from pymotifs.cli import introspect as intro

from pymotifs.utils import flatten
from pymotifs.utils import toposort as topo


class Dispatcher(object):
    """A class which loads and runs stages for the pipeline. This manages
    finding and loading them. It will determine the dependecies for the stage
    and run them as well if desired.
    """

    def __init__(self, name, *args, **kwargs):
        """Create a new dispatcher.

        Parameters
        ----------
        name : str
            String name of the stage to run. Should be like 'update', which
            will load pymotifs.update.
        *args : obj
            Arguments to build the stage with
        skip_dependencies : bool, optional
            Flag to indicate if all dependencies should be skipped. Defaults to
            False.
        exclude : list
            A list, set or tuple of stage names to exclude. This will also
            exclude all dependencies of the stage if they are only used for the
            stage.
        """

        self.name = name
        self._args = args
        self.skip_dependencies = kwargs.get('skip_dependencies')
        self.exclude = set(kwargs.get('exclude', []) or [])
        self.logger = logging.getLogger(__name__)

    def to_exclude(self, *names, **kwargs):
        """Compute a set of stages to exclude. This will load all stages in the
        exclude property to do so. If one entry there is a StageContainer then
        all stages in it will be excluded.
        """

        skip_dependencies = kwargs.get('skip_dependencies', False)
        if skip_dependencies:
            return True

        exclude = set()
        for name in names:
            klass = intro.get_loader(name)
            stage = klass(*self._args)
            exclude.add(name)
            if isinstance(stage, core.StageContainer):
                exclude.update(s.name for s in stage.expand())
        return exclude

    def dependencies(self, stages):
        """Compute the dependency graph for the given stages.

        Parameters
        ----------
        stages : list
            The list of stages to process.

        Returns
        --------
        deps :dict
            A dictonary with the dependency graph.
        """

        def process_stage(stage):
            if issubclass(stage, core.StageContainer):
                return [process_stage(s) for s in stage.stages]
            return [stage]

        stack = []
        for stage in stages:
            if issubclass(stage, core.StageContainer):
                stack.extend(stage.stages)
            else:
                stack.append(stage)

        deps = {}
        while stack:
            current = stack.pop()
            if current in deps:
                next

            if issubclass(current, core.StageContainer):
                stack.extend(current.stages)
            else:
                current_deps = [process_stage(s) for s in current.dependencies]
                current_deps = flatten(current_deps)
                current_deps = set(current_deps)
                deps[current] = current_deps
                stack.extend(current_deps)

        return deps

    def levels(self, dependencies, exclude, allowed):
        """Compute all the stages at each level to run. This will filter the
        stages to only those requested by exclude and allow as well as order
        them. The ordering makes it easier to test and compare if all stages
        are being loaded correctly.
        """

        for level in topo.levels(dependencies):
            current = []
            stages = []
            for k in level:
                try:
                    stages.append(k(*self._args))
                except:
                    self.logger.error("Failed building stage %s", k)
                    raise

            # self.logger.info("Level %s" % level)

            for stage in sorted(stages, key=lambda c: c.name):

                # self.logger.info("Stage %s" % (stage))

                name = stage.name
                if exclude is True or name in exclude:
                    is_allowed = name in allowed or stage in allowed
                    if is_allowed:
                        current.append(stage)
                else:
                    current.append(stage)

            if current:
                yield current

    def flatten(self, dependencies, exclude, allowed):
        """Sort the dependencies into a list of stages to run while respecting
        the exclude and allowed flags. In addition, the stages will be built if
        requested.

        :param dict dependencies: The dependency graph from dependencies.
        :param list exclude: A list of strages to exclude.
        :param bool/list allowed: Either True to indicate all stages are
        allowed or a list of stages which are allowed.
        :param bool build: A flag to indicate if the stage should be built or
        not.
        :returns: A list of the stages as specified by the dependencies and
        exclude/allowed flags.
        """

        levels = self.levels(dependencies, exclude, allowed)
        stages = list(it.chain.from_iterable(levels))
        if not stages:
            raise core.InvalidState("No stages to run")
        return stages

    def stages(self, name):
        """Determine all stages to run and in what order for the given stage
        name. If dependencies is set to True then this will go through all
        dependecies of the given stage and place them in a tree, as well as all
        of their dependecies and so forth. The stages will be sorted
        topologically and then returned in that order.

        If dependecies is False, then a list of one element, the specified
        stage will be returned.

        :param str name: The name of the stage to run.
        :returns: A list of the stages to run.
        """

        allowed = set()
        exclude = self.to_exclude(*self.exclude)
        klass = intro.get_loader(name)
        stage = klass(*self._args)

        if self.skip_dependencies:
            allowed.add(stage.name)
            exclude = True

        deps = None
        if isinstance(stage, core.StageContainer):
            allowed.update(s.name for s in stage.expand())
            allowed.discard(stage.name)
            if exclude is not True:
                allowed.difference_update(exclude)
            if not self.skip_dependencies:
                exclude.add(stage.name)
            deps = self.dependencies(stage.stages)
        else:
            deps = self.dependencies([klass])

        return self.flatten(deps, exclude, allowed)

    def __call__(self, entries, **kwargs):
        """Call the specified stages using the given entries as input. This
        will determine what stages to run using the name property and then run
        them in the correct order.

        :param list entries: The entries to use as input.
        :kwargs: Keyword arguments to pass to each stage.
        """

        stages = self.stages(self.name)
        self.logger.info('Running stages: %s',
                         ', '.join(s.name for s in stages))

        for stage in stages:
            try:
                self.logger.info("Running stage: %s", stage.name)
                stage(entries, **kwargs)
            except Exception as err:
                self.logger.error("Uncaught exception with stage: %s",
                                  self.name)
                raise err
            if 'ml_release_id' in kwargs['manual']:
                del kwargs['manual']['ml_release_id']

        self.logger.info("Finished pipeline")
        print("Finished pipeline")
