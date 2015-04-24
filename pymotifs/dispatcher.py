import logging
import inspect

from pymotifs.core import Stage

from pymotifs.utils.toposort import toposort


class Dispatcher(object):
    """A class which loads and runs stages for the pipeline. This manages
    finding and loading them. It will determine the dependecies for the stage
    and run them as well if desired.
    """

    def __init__(self, name, *args, **kwargs):
        self.name = name
        self._args = args
        self.skip_dependencies = kwargs.get('skip_dependencies')
        self.logger = logging.getLogger(__name__)

    def get_stage(self, name):
        import_name = 'pymotifs.' + name
        from_list = import_name.split('.')
        module = __import__(import_name, fromlist=from_list)
        pairs = inspect.getmembers(module, self.is_loader(import_name))
        if not pairs:
            raise ImportError("Could not find a class to use from: %s" % name)
        return pairs[0][1]

    def stages(self, name):
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

        deps = {stage: stage.dependencies}
        stack = list(stage.dependencies)
        while stack:
            current = stack.pop()
            if current not in deps:
                deps[current] = current.dependencies
                stack.extend(current.dependencies)

        stages = list(toposort(deps))
        if stage not in stages:
            self.logger.warning("Likely there is an issue in toposort")
            stages.append(stage)
        return stages

    def is_loader(self, name):
        def checker(obj):
            return inspect.isclass(obj) and obj.__module__ == name and \
                issubclass(obj, Stage)
        return checker

    def __call__(self, entries, **kwargs):
        try:
            for stage in self.stages(self.name):
                self.logger.info("Running stage: %s", stage)
                obj = stage(*self._args)
                obj(entries, **kwargs)
        except Exception as err:
            self.logger.error("Uncaught exception with stage: %s", self.name)
            self.logger.error("Message: %s" % str(err))
            self.logger.exception(err)
            raise err
