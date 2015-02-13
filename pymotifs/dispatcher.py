import logging
import inspect

from pymotifs.core import Stage


class Dispatcher(object):
    """A class which loads and runs stages for the pipeline. This manages
    finding and loading them. Future work will probably extend this to examine
    and dependecies between the stages and run them in the correct order.
    """

    def __init__(self, name, *args):
        self.name = name
        self._args = args
        self.logger = logging.getLogger(__name__)

    def get_stage(self, name):
        import_name = 'pymotifs.' + name
        from_list = import_name.split('.')
        module = __import__(import_name, fromlist=from_list)
        pairs = inspect.getmembers(module, self.is_loader(import_name))
        if not pairs:
            raise ImportError("Could not find a class to use from: %s" % name)
        return pairs[0][1]

    def is_loader(self, name):
        def checker(obj):
            return inspect.isclass(obj) and obj.__module__ == name and \
                issubclass(obj, Stage)
        return checker

    def __call__(self, entries, **kwargs):
        try:
            stage = self.get_stage(self.name)
            obj = stage(*self._args)
            obj(entries, **kwargs)
        except Exception as err:
            self.logger.error("Uncaught exception with stage: %s", self.name)
            self.logger.error("Message: %s" % str(err))
            self.logger.exception(err)
            raise err
