# Makefile for Arnold aaOcean shader

# this makefile requires openexr, ilmbase and zlib
# please visit the following URL to see how to compile these
# https://bitbucket.org/amaanakram/aaocean/wiki/Home

# $(DEPLOC) is where you have installed 'dependencies' like openexr, ilmbase, and zlib
# the following folder structure is expected under $(DEPLOC)
# $(DEPLOC)/include
# $(DEPLOC)/lib
DEPLOC		:= ../../../../../../sdks/deploy

# this makefile expects the following folder structure under $(ARNOLD) folder
# $(ARNOLD)/include
# $(ARNOLD)/lib
ARNOLD	 	:= $(DEPLOC)/arnold/Arnold-4.2.0.5-linux

####### NO USER-INPUT REQUIRED BEYOND THIS LINE ######

LIBSLOC		:= $(DEPLOC)/lib
GXXLOC	 	:= 
TARGETLOC	:= 

INCLUDES := -I$(ARNOLD)/include -I../../../../externals/aaOcean/src -I../../../../externals/helpers  -I$(DEPLOC)/include -I$(DEPLOC)/include/OpenEXR
LDFLAGS  := -fopenmp -L$(ARNOLD)/bin 
LIBS     := $(LIBSLOC)/libIlmImf.a $(LIBSLOC)/libImath.a $(LIBSLOC)/libHalf.a $(LIBSLOC)/libIexMath.a $(LIBSLOC)/libIlmThread.a $(LIBSLOC)/libIex.a $(LIBSLOC)/libz.a -Wl,-rpath -Wl,$(LIBSLOC) -lai

CXXFLAGS := -w -fPIC -fopenmp -fvisibility=hidden -pthread -DWRITE_OPENEXR -ffast-math -fno-strict-aliasing -msse2 -O3

CXX	 := $(GXXLOC)g++
SOURCES  := shaders/aaOceanArnoldShader.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := $(TARGETLOC)aaOceanArnold.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS) -L$(LIBSLOC) $(LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
