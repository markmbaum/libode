#configuration file, which specifies compiler, flags, paths, etc.
include config.mk

#-------------------------------------------------------------------------------
#stuff to compile

#objects to be built
other=ode_misc.o
rkbase=OdeBase.o \
		OdeBaseRK.o \
		OdeBaseARK.o
rkmeth=OdeEuler.o \
		OdeMidpnt.o \
		OdeTrapz.o \
		OdeRKF21.o \
		OdeSsp3.o \
		OdeRKF32.o \
		OdeRK4.o \
		OdeRK43.o \
		OdeRKCK.o \
		OdeDoPri54.o \
		OdeButcher6.o \
		OdeDoPri87.o

#executables to be built for testing (if `make tests` is called)
execs=ode_test_adapt.exe \
		ode_test_fixed.exe \
		ode_test_snaps.exe \
		ode_test_conv.exe \
		ode_test_work.exe \
		ode_test_swarmalator.exe

#-------------------------------------------------------------------------------
#prepend directories to names

oth=$(patsubst %,$(diro)/%,$(other))
rkb=$(patsubst %,$(diro)/%,$(rkbase))
rkm=$(patsubst %,$(diro)/%,$(rkmeth))
exe=$(patsubst %,$(dirb)/%,$(execs))

#-------------------------------------------------------------------------------

#automatic target for everything except the tests
all: $(dirb)/libode.a

#target for everything including tests
tests: all $(exe)

#-------------------------------------------------------------------------------
#compile objects

$(oth): $(diro)/%.o: $(dirs)/%.cpp $(dirs)/%.hpp
	$(cxx) $(cflags) -o $@ -c $<

$(rkb): $(diro)/%.o: $(dirs)/%.cpp $(dirs)/%.hpp $(oth)
	$(cxx) $(cflags) -o $@ -c $<

$(rkm): $(diro)/%.o: $(dirs)/%.cpp $(dirs)/%.hpp $(oth) $(rkb)
	$(cxx) $(cflags) -o $@ -c $<

#-------------------------------------------------------------------------------
#create library of solvers from the compiled objects

$(dirb)/libode.a: $(rkb) $(rkm) $(oth)
	ar r $(dirb)/libode.a $(rkm) $(rkb) $(oth)

#-------------------------------------------------------------------------------
#compile test programs

$(exe): $(dirb)/%.exe: $(dirt)/%.cpp $(rkm) $(dirt)/ode_explicit_test_systems.hpp
	$(cxx) $(cflags) -o $@ $< -I $(dirs) -L$(dirb) -lode

#-------------------------------------------------------------------------------
#other things

#target to remove all the built files
clean:
	rm  $(diro)/*.o $(dirb)/*.exe $(dirb)/*.a

.PHONY: clean
