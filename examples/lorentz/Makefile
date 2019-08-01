#main libode directory
dtop=../..

#get configuration from the top directory config file
include $(dtop)/config.mk

#directory with header files
I=-I$(dtop)/$(ds)

all: libodemake lorentz.exe

#rule for jumping to the libode makefile
libodemake:
	$(MAKE) -C $(dtop)

#rule for compiling the integrator
lorentz.exe: lorentz.cc lorentz.h libodemake
	$(cxx) $(flags) -o $@ $< $(I) -L$(dtop)/$(db) -lode
