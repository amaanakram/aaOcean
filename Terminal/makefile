# makefile for aaOcean Terminal/Shell application

# this makefile requires openexr, ilmbase and zlib
# please visit the following URL to see how to compile these
# https://bitbucket.org/amaanakram/aaocean/wiki/Home

# $(DEPLOC) is where you have installed 'dependencies' like openexr, ilmbase, and zlib
# the following folder structure is expected under $(DEPLOC)
# $(DEPLOC)/include
# $(DEPLOC)/lib
DEPLOC      := ../../../sdks/deploy

####### NO USER-INPUT REQUIRED BEYOND THIS LINE ######

AAOCEANLOC    := ../../aaocean

# the compiler: gcc for C program, define as g++ for C++
CXX     := $(GXXLOC)g++

# compiler flags:
CXXFLAGS := -w -ffast-math -fPIC -fopenmp -fvisibility=hidden -pthread -fno-strict-aliasing -msse2 -O3

# includes
INCLUDES := -I$(AAOCEANLOC)/externals/aaOcean/src -I$(AAOCEANLOC)/externals/helpers -I$(DEPLOC)/include -I$(DEPLOC)/include/OpenEXR

# libs where exr libs and zlib are installed:
LIBSLOC     := $(DEPLOC)/lib

# libs to link against 
LIBS     := $(LIBSLOC)/libIlmImf.a $(LIBSLOC)/libImath.a $(LIBSLOC)/libHalf.a $(LIBSLOC)/libIexMath.a $(LIBSLOC)/libIlmThread.a $(LIBSLOC)/libIex.a $(LIBSLOC)/libz.a -Wl,-rpath -Wl,$(LIBSLOC)

# this should not be changed
LDFLAGS  := -fPIC -fopenmp -pthread -fvisibility=hidden 

# the build target executable:
TARGET = aaOceanTerminal

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(TARGET) $(TARGET).cpp $(LDFLAGS) -L$(LIBSLOC) $(LIBS)

clean:
	$(RM) $(TARGET)
