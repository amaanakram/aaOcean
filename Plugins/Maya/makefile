# makefile for maya aaocean deformer

# https://bitbucket.org/amaanakram/aaocean/wiki/Home

# This makefile is for gcc-4.2.x, which is ABI-compatible with gcc-4.1.2
# makefile will also work for subsequent gcc versions, but please check
# ABI compatibility with your maya version

# this makefile expects the following folder structure under $(MAYA) folder
# $(MAYA)/include
# $(MAYA)/lib

# $(DEPLOC) is where you have installed 'dependencies' like openexr, ilmbase, and zlib
# the following folder structure is expected under $(DEPLOC)
# $(DEPLOC)/include
# $(DEPLOC)/lib

DEPLOC		:= ../../../../sdks/deploy
MAYA		:= $(DEPLOC)/maya/maya2014

GXXLOC		:= #path to gcc here
TARGETLOC	:= 

####### NO USER-INPUT REQUIRED BEYOND THIS LINE ######

MAYALIB		:= $(MAYA)/lib
LIBSLOC		:= $(DEPLOC)/lib

INCLUDES	:= -I$(MAYA)/include -I$(DEPLOC)/include -I../../externals/aaOcean/src -I../../externals/helpers
LDFLAGS		:= -fopenmp
LIBS		:= -L$(MAYALIB) -lOpenMaya -lOpenMayaAnim
CXXFLAGS	:= -w -O3 -fPIC -fvisibility=hidden -DLINUX -D_BOOL -DI386 -ffast-math -DREQUIRE_IOSTREAM -fno-strict-aliasing -msse2 
CXX		:= $(GXXLOC)g++
SOURCES		:= aaOceanDeformer.cpp
OBJECTS		:= $(SOURCES:.cpp=.o)
TARGET		:= $(TARGETLOC)aaOceanDeformer.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS) -L$(LIBSLOC) $(LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
