#TEMPLATE CONFIG FILE

#-------------------------------------------------------------------------------
#choose a compiler

#c++ compiler
cxx=g++ # <--- GNU C++ compiler
#cxx=icpc # <--- Intel C++ compiler

#-------------------------------------------------------------------------------
#choose compilation flags

#compilation flags
cflags=-std=c++11 -Wall -Wextra -pedantic -O3 # <--- GNU compiler flags
#cflags=-std=c++11 -Wall -O3 # <--- Intel compiler flags

#-------------------------------------------------------------------------------
#provide directory names (these probably shouldn't be changed)

#source code directory
dirs=src
#test code directory
dirt=test/cpp
#compiled object directory
diro=obj
#compiled executable directory
dirb=bin
