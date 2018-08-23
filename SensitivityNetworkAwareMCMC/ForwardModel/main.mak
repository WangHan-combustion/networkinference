
# This Makefile builds a C++ application that uses Cantera.  By
# default, the main program file is 'main.cpp,' which prints out some
# properties of a reacting gas mixture. 

# To build program 'main', simply type 'make', or 'make -f <this
# file>' if this file is named something other than 'Makefile.'  

# Once you have verified that the demo runs, edit this file to replace
# object file 'main.o' with your own object file or files. 


#------------------------  edit this block ---------------------------------

# the name of the executable program to be created
PROG_NAME = main

# the object files to be linked together. 
OBJS = main.o SOFCModel.o SOFCAnode.o SOFCAnodeIntegrator.o gauss_legendre.o structDef.o  fileRead.o fileWrite.o gaussPatterson.o tools.o

# additional flags to be passed to the linker. If your program
# requires other external libraries, put them here
LINK_OPTIONS =   -L/usr/local/lib 

#---------------------------------------------------------------------------
# You probably don't need to edit anything below.

# the C++ compiler
CXX = gcc 

# C++ compile flags
CXXFLAGS = -O3 -Wall   

# external libraries
EXT_LIBS = -luser -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo -lctnumerics -lctmath -ltpx -lctspectra -lconverters -lctbase -lsundials_cvodes -lsundials_nvecserial -lsundials_ida -lpthread -lctf2c -lctcxx -lgfortran -lgsl -lgslcblas -llapack -lm -lblas 

# Ending C++ linking libraries
LCXX_END_LIBS = -lctf2c -lm -lstdc++ 

# the directory where the Cantera libraries are located
CANTERA_LIBDIR=/Applications/Cantera/lib 

# the directory where Cantera include files may be found.
CANTERA_INCDIR=/Applications/Cantera/include

# flags passed to the C++ compiler/linker for the linking step
LCXXFLAGS = -L$(CANTERA_LIBDIR)  -L/opt/intel/Compiler/11.1/073/mkl/lib/em64t -L/usr/local/lib -O3 -Wall  

# C pre-processor flags
CPPFLAGS += -I$(CANTERA_INCDIR)

PROGRAM = $(PROG_NAME)$(EXE_EXT) 

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CXX) -o $(PROGRAM) $(OBJS) $(LCXXFLAGS)\
        $(CANTERA_LIBS) $(LINK_OPTIONS) $(EXT_LIBS) \
          $(LCXX_END_LIBS) 

clean:
	rm -f *.d*
	rm -f $(OBJS)
	rm -f $(PROGRAM)

TAGS: 
	etags *.h *.cpp

.PHONY: clean all

# Scan for dependencies if not cleaning
ifneq ($(MAKECMDGOALS), clean)
  include $(subst .o,.d,$(notdir $(OBJS)))
endif

# Depdendency file generation - note that these files are all an extension of 
# the same command using \ (so only one running shell) to ensure that the 
# process ID $$ stays the same.
%.d: %.cpp
	$(call dependency-generation,$(CXX))

%.d: %.c
	$(call dependency-generation,$(CC))

# Function argument is compiler name
define dependency-generation
	$1 $(CPPFLAGS) -MM -MF $@.$$$$ -MT $(TARGET_DIR)/$*.o $<; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@;       \
	rm -f $@.$$$$
endef
