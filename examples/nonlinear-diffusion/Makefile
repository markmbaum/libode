#main libode directory
dtop=../..

#get configuration from the top directory config file
include $(dtop)/config.mk

#directory with header files
I=-I$(dtop)/$(ds)

all: libodemake nonlinear_diffusion.exe

#rule for jumping to the libode makefile
libodemake:
	$(MAKE) -C $(dtop)

#rule for compiling the integrator
nonlinear_diffusion.exe: nonlinear_diffusion.cc nonlinear_diffusion.h libodemake
	$(cxx) $(flags) -o $@ $< $(I) -L$(dtop)/$(db) -lode
