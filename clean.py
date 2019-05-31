#short script to remove compiled objects, executables, grid files, and model output

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
