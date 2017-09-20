"""From:
http://scottlobdell.me/2015/08/using-decorators-python-automatic-registration/
"""
REGISTERED_CLASSES = {}

def register_class(cls):
    REGISTERED_CLASSES[lower(cls.__name__)] = cls
    return cls
