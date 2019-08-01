#main libode directory
dtop=../..

#get configuration from the top directory config file
include $(dtop)/config.mk

#directory with libode header files
I=-I$(dtop)/$(ds)

all: libodemake swarmalator.exe

#rule for jumping to the libode makefile
libodemake:
	$(MAKE) -C $(dtop)

#rule for compiling the integrator
swarmalator.exe: swarmalator.cc swarmalator.h libodemake
	$(cxx) $(flags) -o $@ $< $(I) -L$(dtop)/$(db) -lode
