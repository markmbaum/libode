#name of example
name=swarmalator
#main libode directory
dtop=../..
#get configuration from the top directory config file
include $(dtop)/config.mk

#directory with header files
I=-I$(dtop)/include/ode

#directory with libode library
L=-L$(dtop)/bin

all: libodemake $(name).exe

#rule for jumping to the libode makefile
libodemake:
	$(MAKE) -C $(dtop)

#rule for compiling the integrator
$(name).exe: $(name).cc $(name).h libodemake
	$(CXX) $(CFLAGS) -o $@ $< $(I) $(L) -lode
