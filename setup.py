# script.py
import argparse
from subprocess import call

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('cmd', type=str, help='setup command')
    
    args = parser.parse_args()
    cmd = args.cmd
    
    # magellan
    call(["python", "magellan/setup.py", cmd])
    
    
