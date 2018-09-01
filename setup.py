# setup wrapper for unicorn
from subprocess import call
import sys
import glob
import os


if __name__ == '__main__':

    args = [arg for arg in sys.argv[1:]]
   
    for module in glob.glob('*'):
        if os.path.isdir(module):
            if ('docs' in module)\
               or ('build' in module)\
               or ('dist' in module)\
               or ('egg' in module):
                continue

            print '-------> trying to operate on module: ', module

            try:
                call(
                    ["python", "setup.py"] + args,
                    cwd=module)
            except Exception:
                print 'cannot operate on module %s' % module

    
    
