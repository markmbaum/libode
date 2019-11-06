#TEMPLATE CONFIG FILE

#-------------------------------------------------------------------------------
#choose a compiler

#c++ compiler
cxx=g++ # <--- GNU C++ compiler
#cxx=icpc # <--- Intel C++ compiler

#-------------------------------------------------------------------------------
#choose compilation flags

#compilation flags
flags=-Wall -Wextra -pedantic -O3 # <--- GNU compiler flags
#flags=-Wall -O3 # <--- Intel compiler flags

#-------------------------------------------------------------------------------
#local directory names (don't change these!)

#source code directory
ds=src
#test code directory
dt=test
#compiled object directory
do=obj
#compiled executable directory
db=bin
