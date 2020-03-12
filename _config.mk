#TEMPLATE CONFIG FILE
#Copy or rename this file "config.mk" and edit the variables below as needed

#-------------------------------------------------------------------------------
#choose a compiler

cxx=g++ # <--- GNU C++ compiler
#cxx=icpc # <--- Intel C++ compiler

#-------------------------------------------------------------------------------
#choose compilation flags

flags=-Wall -Wextra -pedantic -O3 # <--- GNU compiler flags
#flags=-w2 -O3 # <--- Intel compiler flags

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#local directory names, DON'T CHANGE THESE

#source code directory
ds=src
#test code directory
dt=tests
#compiled object directory
do=obj
#compiled executable directory
db=bin
