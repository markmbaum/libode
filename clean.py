"""
This is a short script to remove compiled objects and executables.
It removes all non-hidden files in the 'obj', 'bin', and 'out' directories.
To specify which of these folders to target, use folder names at the command line.
"""

import os
import sys
from os.path import join

if len(sys.argv) > 1:
    dirs = [str(arg) for arg in sys.argv[1:]]
else:
    dirs = ['obj', 'bin', 'out']

for d in dirs:
    for fn in os.listdir(d):
        if fn[0] != '.':
            os.remove(join(d, fn))
    print('directory "%s" cleaned' % d)
