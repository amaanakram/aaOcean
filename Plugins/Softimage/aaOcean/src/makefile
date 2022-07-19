# makefile for softimage aaocean deformer

# this makefile requires fftw library
# please visit the following URL to see how to compile fftw
# https://bitbucket.org/amaanakram/aaocean/wiki/Home

# this makefile expects the following folder structure under $(SOFTIMAGE) folder
# $(SOFTIMAGE)/include
# $(SOFTIMAGE)/lib

# $(DEPLOC) is where you have installed 'dependencies' like openexr, ilmbase, and zlib
# the following folder structure is expected under $(DEPLOC)
# $(DEPLOC)/include
# $(DEPLOC)/lib

DEPLOC		:= ../../../../../../sdks/deploy
SOFTIMAGE	:= $(DEPLOC)/softimage/softimage2013

####### NO USER-INPUT REQUIRED BEYOND THIS LINE ######

LIBSLOC		:= $(DEPLOC)/lib
GXXLOC		:= 
TARGETLOC	:= 

INCLUDES 	:= -I$(SOFTIMAGE)/include -I$(DEPLOC)/include -I../../../../externals/aaOcean/src -I../../../../externals/helpers
LDFLAGS		:= -fopenmp
LIBS		:=
CXXFLAGS	:= -w -O3 -fPIC -ffast-math -fno-strict-aliasing -msse2
CXX		:= $(GXXLOC)g++
SOURCES  	:= icenode.cpp
OBJECTS  	:= $(SOURCES:.cpp=.o)
TARGET   	:= $(TARGETLOC)aaOceanICE.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS) -L$(LIBSLOC) $(LIBS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
