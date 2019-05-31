import sys
from os import rename, remove
from os.path import join, isfile

sys.path.append(join('..', '..', 'odpydoc'))
from odpydoc import doc

sys.path.append(join('..'))
doc('pylibode', outdir=join('..', 'docs'))
