import logging
import inspect

from pymotifs.core import Stage
from pymotifs.core import InvalidState
from pymotifs.core import StageContainer
from pymotifs.cli import introspect as intro

from pymotifs.utils.toposort import toposort


class Dispatcher(object):
    """A class which loads and runs stages for the pipeline. This manages
    finding and loading them. It will determine the dependecies for the stage
    and run them as well if desired.
    """

    def __init__(self, name, *args, **kwargs):
        """Create a new dispatcher.

        :name: String name of the stage to run. Should be like 'update', which
        will load pymotifs.update.
        :*args: Arguments to build the stage with
        :skip_dependencies: Flag to indicate if all dependencies should be
        skipped.
        :exclude: A list, set or tuple of stage names to exclude. This will
        also exclude all dependencies of the stage if they are only used for
        the stage.
        """

        self.name = name
        self._args = args
        self.skip_dependencies = kwargs.get('skip_dependencies')
        self.exclude = set(kwargs.get('exclude', []) or [])
        self.logger = logging.getLogger(__name__)

    def to_exclude(self):
        """Compute a set of stages to exclude. This will load all stages in the
        exclude property to do so. If one entry there is a StageContainer then
        all stages in it will be excluded.
        """
        exclude = set()
        for name in self.exclude:
            stage = intro.get_loader(name)
            exclude.add(stage)
            if issubclass(stage, StageContainer):
                exclude.update(stage.stages)
        return exclude

    def __builder__(self, build):
        """Create a function to build the stages. If build is True it will
        construct the stages otherwise it will return the given stage. If given
        a StageContainer it will return a list of all stages in the
        StageContainer.

        :param bool build: Flag to indicate if the stage should be built.
        :returns: A function for building stages.
        """

        def fn(stage):
            if issubclass(stage, StageContainer):
                stages = []
                for s in stage.stages:
                    stages.extend(fn(s))
                return stages
            if build:
                return [stage(*self._args)]
            return [stage]
        return fn

    def stages(self, name, build=False):
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

        fn = self.__builder__(build)
        stage = intro.get_loader(name)
        if self.skip_dependencies:
            return fn(stage)

        exclude = self.to_exclude()
        deps = {stage: stage.dependencies}
        stack = list(stage.dependencies)

        if issubclass(stage, StageContainer):
            deps = {}
            exclude.add(stage)
            stack = list(stage.stages)

        while stack:
            current = stack.pop()
            if current in exclude or current in deps:
                next

            to_add = getattr(current, 'dependencies', [])
            if issubclass(current, StageContainer):
                exclude.add(current)
                to_add = getattr(current, 'stages')

            deps[current] = to_add
            stack.extend(to_add)

        stages = []
        for s in toposort(deps, by=lambda s: s.__name__):
            if s not in exclude:
                stages.extend(fn(s))

        if not stages:
            raise InvalidState("No stages to run")

        return stages

    def __call__(self, entries, **kwargs):
        stages = self.stages(self.name, build=True)
        names = [s.name for s in stages]
        self.logger.debug('Running stages: %s', ', '.join(names))
        for stage in stages:
            try:
                self.logger.info("Running stage: %s", stage.name)
                stage(entries, **kwargs)
            except Exception as err:
                self.logger.error("Uncaught exception with stage: %s",
                                  self.name)
                self.logger.error("Message: %s" % str(err))
                self.logger.exception(err)
                raise err
