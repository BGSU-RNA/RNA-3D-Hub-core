import logging
import inspect

from pymotifs.core import Stage
from pymotifs.core import InvalidState
from pymotifs.core import StageContainer

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
        self.exclude = set(kwargs.get('exclude', []))
        self.logger = logging.getLogger(__name__)

    def get_stage(self, name):
        import_name = 'pymotifs.' + name
        from_list = import_name.split('.')
        module = __import__(import_name, fromlist=from_list)
        pairs = inspect.getmembers(module, self.is_loader(import_name))
        if not pairs:
            raise ImportError("Could not find a class to use from: %s" % name)
        return pairs[0][1]

    def to_exclude(self):
        """Compute a set of stages to exclude. This will load all stages in the
        exclude property to do so. If one entry there is a StageContainer then
        all stages in it will be excluded.
        """
        exclude = set()
        for name in self.exclude:
            stage = self.get_stage(name)
            exclude.add(stage)
            if issubclass(stage, StageContainer):
                exclude.update(stage.stages)
        return exclude

    def stages(self, name, build=False):
        """Determine all stages to run and in what order for the given stage
        name. If dependencies is set to True then this will go through all
        dependecies of the given stage and place them in a tree, as well as all
        of their dependecies and so forth. The stages will be sorted
        topologically and then returned in that order.

        If dependecies is False, then a list of one element, the specified
        stage will be returned.

        :param str name: The name of the stage to run.
        """

        stage = self.get_stage(name)
        if self.skip_dependencies:
            return [stage]

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

        fn = lambda s: s
        if build:
            fn = lambda s: s(*self._args)

        by = lambda s: s.__name__
        stages = [fn(s) for s in toposort(deps, by=by) if s not in exclude]

        if stage not in stages and stage not in exclude:
            self.logger.warning("Likely there is an issue in toposort")
            stages.append(fn(stage))

        if not stages:
            raise InvalidState("No stages to run")

        return stages

    def is_loader(self, name):
        def checker(obj):
            return inspect.isclass(obj) and obj.__module__ == name and \
                issubclass(obj, Stage)
        return checker

    def __call__(self, entries, **kwargs):
            for stage in self.stages(self.name, build=True):
                try:
                    self.logger.info("Running stage: %s", stage)
                    stage(entries, **kwargs)
                except Exception as err:
                    self.logger.error("Uncaught exception with stage: %s",
                                      self.name)
                    self.logger.error("Message: %s" % str(err))
                    self.logger.exception(err)
                    raise err
