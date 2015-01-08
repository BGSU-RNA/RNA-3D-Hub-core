import logging
import inspect


class Dispatcher(object):
    """A class which loads and runs stages for the pipeline. This manages
    finding and loading them. Future work will probably extend this to examine
    and dependecies between the stages and run them in the correct order.
    """

    def __init__(self, name, *args):
        self.stage = self.get_stage(name)
        self._args = args
        self.logger = logging.getLogger(__name__)

    def get_stage(self, name):
        fromlist = ['pymotifs']
        if '.' in name:
            fromlist.extend(name.split('.')[1:])
        module = __import__(name, fromlist=fromlist)
        pairs = inspect.getmembers(module, self.is_loader(name))
        if not pairs:
            raise ValueError("Could not find a class to use from: %s" % name)
        return pairs[0][1]

    def is_loader(self, name):
        def checker(obj):
            return inspect.isclass(obj) and obj.__module__ == name
        return checker

    def __call__(self, entries, **kwargs):
        try:
            obj = self.stage(*self._args)
            obj(entries)
        except Exception as err:
            self.logger.error("Uncaught exception with stage: %s", self.stage)
            self.logger.exception(err)
