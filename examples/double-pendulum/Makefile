#main libode directory
dtop=../..

#get configuration from the top directory config file
include $(dtop)/config.mk

#directory with header files
I=-I$(dtop)/$(ds)

all: libodemake double_pendulum.exe

#rule for jumping to the libode makefile
libodemake:
	$(MAKE) -C $(dtop)

#rule for compiling the integrator
double_pendulum.exe: double_pendulum.cc double_pendulum.h libodemake
	$(cxx) $(flags) -o $@ $< $(I) -L$(dtop)/$(db) -lode
