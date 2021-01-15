#you probably shouldn't need to edit anything in here

-include config.mk

#-------------------------------------------------------------------------------
#local directory names

#source code directory
ds=src
#test code directory
dt=tests
#compiled object directory
do=obj
#compiled executable directory
db=bin

#-------------------------------------------------------------------------------
#stuff to compile

#modules with support functions or classes
auxil=  ode_io.o \
		ode_util.o \
		ode_linalg.o

#base classes
base=   ode_base.o \
		ode_adaptive.o \
		ode_embedded.o \
		ode_rk.o \
		ode_erk.o \
		ode_irk.o \
		ode_rosenbrock.o \
		ode_newton.o \

#solver classes
solver= ode_euler.o \
		ode_trapz.o \
		ode_ssp_3.o \
		ode_rkf_32.o \
		ode_rk_4.o \
		ode_rk_43.o \
		ode_dopri_54.o \
		ode_rkck.o \
		ode_vern_65.o \
		ode_vern_76.o \
		ode_dopri_87.o \
		ode_vern_98.o \
		ode_grk4a.o \
		ode_row6a.o \
		ode_backward_euler.o \
		ode_gauss_6.o \
		ode_lobatto_iiic_6.o \
		ode_radau_iia_5.o \
		ode_geng_5.o \
		ode_sdirk_43.o

#executables to be built for testing when `make tests` is used
execs=  test_fixed.exe \
		test_all.exe \
		test_snaps.exe \
		test_conv.exe \
		test_adapt.exe \
		test_stiff.exe \
		test_work.exe \
		test_linalg.exe \
		test_newton.exe

#-------------------------------------------------------------------------------
#prepend appropriate directories to names

aux=$(patsubst %,$(do)/%,$(auxil))
bas=$(patsubst %,$(do)/%,$(base))
sol=$(patsubst %,$(do)/%,$(solver))
exe=$(patsubst %,$(db)/%,$(execs))

#-------------------------------------------------------------------------------
#make some short names for the base classes

B=ode_base
BA=ode_adaptive
BE=ode_embedded
BRK=ode_rk
BERK=ode_erk
BIRK=ode_irk
BR=ode_rosenbrock
N=ode_newton
NB=ode_newton_bridge

#-------------------------------------------------------------------------------

#automatic target for everything except the tests
all: $(db)/libode.a

#target for everything including tests
tests: all $(exe)

#-------------------------------------------------------------------------------
#compile objects

#auxils, the support functions and stuff don't have complicated dependencies
$(aux): $(do)/%.o: $(ds)/%.cc $(ds)/%.h
	$(CXX) $(CFLAGS) -o $@ -c $<

#OdeNewton
$(do)/$(N).o: $(ds)/$(N).cc $(ds)/$(N).h $(ds)/$(NB).h $(do)/ode_linalg.o
	$(CXX) $(CFLAGS) -o $@ -c $< -I$(ds)

#OdeBase
$(do)/$(B).o: $(ds)/$(B).cc $(ds)/$(B).h $(do)/ode_io.o $(do)/ode_util.o
	$(CXX) $(CFLAGS) -o $@ -c $< -I$(ds)
#OdeAdaptive
$(do)/$(BA).o: $(ds)/$(BA).cc $(ds)/$(BA).h $(do)/$(B).o $(do)/ode_io.o $(do)/ode_util.o
	$(CXX) $(CFLAGS) -o $@ -c $< -I$(ds)

#OdeRK
$(do)/$(BRK).o: $(ds)/$(BRK).cc $(ds)/$(BRK).h
	$(CXX) $(CFLAGS) -o $@ -c $< -I$(ds)
#OdeERK
$(do)/$(BERK).o: $(ds)/$(BERK).cc $(ds)/$(BERK).h
	$(CXX) $(CFLAGS) -o $@ -c $< -I$(ds)
#OdeIRK
$(do)/$(BIRK).o: $(ds)/$(BIRK).cc $(ds)/$(BIRK).h
	$(CXX) $(CFLAGS) -o $@ -c $< -I$(ds)
#OdeEmbedded
$(do)/$(BE).o: $(ds)/$(BE).cc $(ds)/$(BE).h $(do)/$(BA).o $(do)/ode_util.o
	$(CXX) $(CFLAGS) -o $@ -c $< -I$(ds)
#OdeRosenbrock
$(do)/$(BR).o: $(ds)/$(BR).cc $(ds)/$(BR).h $(do)/ode_linalg.o
	$(CXX) $(CFLAGS) -o $@ -c $< -I$(ds)

#solver classes can all just depend on the base classes
$(sol): $(do)/%.o: $(ds)/%.cc $(ds)/%.h $(bas)
	$(CXX) $(CFLAGS) -o $@ -c $< -I$(ds)

#-------------------------------------------------------------------------------
#create library of solvers from the compiled objects

$(db)/libode.a: $(bas) $(sol) $(aux)
	ar r $(db)/libode.a $(aux) $(bas) $(sol)

#-------------------------------------------------------------------------------
#compile test programs

$(exe): $(db)/%.exe: $(dt)/%.cc $(sol) $(bas) $(aux) $(dt)/test_systems.h
	$(CXX) $(CFLAGS) -o $@ $< -I$(ds) -L$(db) -lode

#-------------------------------------------------------------------------------
#auxil things

#target to remove all the compiled files
clean:
	-rm $(do)/*.o
	-rm $(db)/*.exe
	-rm $(db)/*.a

.PHONY: clean
