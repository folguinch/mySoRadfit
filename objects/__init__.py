from os.path import basename, dirname, join
from glob import glob

from .register import REGISTERED_CLASSES

# Load all the classes
pwd = dirname(__file__)
for x in glob(join(pwd, '*.py')):
    if not x.startswith('__'):
        __import__(basename(x)[:-3], globals(), locals())

__all__ = ['REGISTERED_CLASSES']
