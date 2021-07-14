# If you are using this Makefile standalone and fastjet-config is not
# in your path, edit this line to specify the full path
FASTJETCONFIG = fastjet-config
PREFIX = $(shell $(FASTJETCONFIG) --prefix)
CXX = g++
CXXFLAGS += -O3 -ffast-math -Wall -std=c++14 -fPIC -DPIC -DNDEBUG -DPIRANHA_USE_PYFJCORE -DEVENTGEOMETRY_USE_PYFJCORE
CXXFLAGS += -IEventGeometry -IEventGeometry/Wasserstein -IEventGeometry/PyFJCore
LDFLAGS += -LEventGeometry -lEventGeometry
install_script = $(SHELL) ./scripts/install-sh
check_script = ./scripts/check.sh

# global contrib-wide Makefile include may override some of the above
# variables (leading "-" means don't give an error if you can't find
# the file)
-include ../.Makefile.inc

#------------------------------------------------------------------------
# things that are specific to this contrib
NAME = Piranha
SRCS = DynamicVoronoiBase.cc DynamicVoronoiCylinder.cc DynamicVoronoiDisk.cc IteratedVoronoiSubtractorBase.cc OptimalTransportSubtractor.cc RecursiveSafeSubtractor.cc
EXAMPLES = example
INSTALLED_HEADERS = DynamicVoronoiBase.hh DynamicVoronoiCylinder.hh DynamicVoronoiDisk.hh IteratedVoronoiSubtractor.hh IteratedVoronoiSubtractorBase.hh OptimalTransportSubtractor.hh Piranha.hh PiranhaUtils.hh RecursiveSafeSubtractor.hh
#------------------------------------------------------------------------

OBJS := $(SRCS:.cc=.o)
EXAMPLES_SRCS  = $(EXAMPLES:=.cc)

install_HEADER  = $(install_script) -c -m 644
install_LIB     = $(install_script) -c -m 644
install_DIR     = $(install_script) -d
install_DATA    = $(install_script) -c -m 644
install_PROGRAM = $(install_script) -c -s
install_SCRIPT  = $(install_script) -c

ifeq "$(shell uname)" "Darwin"
	dynlibopt = -dynamiclib
    dynlibext = dylib
    LDFLAGS += -undefined dynamic_lookup -install_name @rpath/lib$(NAME).dylib -Wl,-rpath,@loader_path/EventGeometry -lgmp
else
    dynlibopt = -shared
    dynlibext = so
    LDFLAGS += -Wl,-rpath,\$$ORIGIN/EventGeometry -lgmp
endif

FJ_CXXFLAGS+= $(shell $(FASTJETCONFIG) --cxxflags)
FJ_LDFLAGS += -lEventGeometry $(shell $(FASTJETCONFIG) --libs) 

.PHONY: clean distclean examples check install all shared install_shared

# http://make.mad-scientist.net/papers/advanced-auto-dependency-generation/#combine
DEPDIR = .deps
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.d
COMPILE.cxx = $(CXX) $(DEPFLAGS) $(CXXFLAGS) -c

%.o : %.cc
%.o : %.cc $(DEPDIR)/%.d | $(DEPDIR)
	$(COMPILE.cxx) $(OUTPUT_OPTION) $<

# compilation of the code (default target)
all: shared

lib$(NAME).a: $(OBJS) 
	ar cru lib$(NAME).a $(OBJS)
	ranlib lib$(NAME).a

shared: lib$(NAME).$(dynlibext)
lib$(NAME).$(dynlibext): $(OBJS)
	$(CXX) $(OBJS) $(dynlibopt) $(LDFLAGS) -o lib$(NAME).$(dynlibext)

# building the examples
examples: $(EXAMPLES)

# pileup comparison
pileup_compare: pileup_compare.o
	$(CXX) $(CXXFLAGS) -o $@ $< -L. -l$(NAME) $(LDFLAGS) -lboost_program_options -lboost_iostreams -lHepPID -lConstituentSubtractor

# the following construct makes it possible to automatically build
# each of the examples listed in $EXAMPLES
$(EXAMPLES): % : %.o all
	$(CXX) $(CXXFLAGS) -o $@ $< -L. -l$(NAME) $(LDFLAGS)

# check that everything went fine
check: examples
	@for prog in $(EXAMPLES); do\
	  $(check_script) $${prog} ../data/single-event.dat || exit 1; \
	done
	@echo "All tests successful"

# cleaning the directory
clean:
	rm -fv *~ *.o *.a *.so *.dylib
	rm -rfv .deps

distclean: clean
	rm -f lib$(NAME).a $(EXAMPLES)

install_shared: shared
	$(install_DIR) $(PREFIX)/lib
	$(install_LIB) lib$(NAME).$(dynlibext) $(PREFIX)/lib

# install things in PREFIX/...
install: all install_shared
	$(install_DIR) $(PREFIX)/include/fastjet/contrib
	for header in $(INSTALLED_HEADERS); do\
	  $(install_HEADER) $$header $(PREFIX)/include/fastjet/contrib/;\
	done
	$(install_DIR) $(PREFIX)/lib
	$(install_LIB) lib$(NAME).a $(PREFIX)/lib

$(DEPDIR): ; @mkdir -p $@

ALL_SRCS = $(SRCS) $(EXAMPLES_SRCS) pileup_compare.cc
DEPFILES := $(ALL_SRCS:%.cc=$(DEPDIR)/%.d)

$(DEPFILES):

include $(wildcard $(DEPFILES))
