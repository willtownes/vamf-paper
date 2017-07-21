"""
Miscellaneous functions used by other scripts and modules
"""
import os,errno

def mkdir_p(path):
    """http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python"""
    try: os.makedirs(path)
    except OSError, exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path): pass
        else: raise exc