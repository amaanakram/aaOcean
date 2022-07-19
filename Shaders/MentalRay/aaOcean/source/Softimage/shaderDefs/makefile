# Makefile for mental ray aaOcean shader's Shader Definition Plugin for Softimage

# this makefile expects the following folder structure under $(SOFTIMAGE) folder
# $(SOFTIMAGE)/include
# $(SOFTIMAGE)/lib
XSISDK 		:= ../../../../../../../../sdks/deploy/softimage/softimage2013
GXXLOC	 	:= 
TARGETLOC	:= 

####### NO USER-INPUT REQUIRED BEYOND THIS LINE ######

INCLUDES := -I$(XSISDK)/include
LDFLAGS  := -Wl,-rpath="\\$$ORIGIN:\\$$ORIGIN/../lib:\\$$ORIGIN/../bin"
CXXFLAGS := -w -O2 -fPIC
CXX      := $(GXXLOC)g++
SOURCES  := shaderDefsPlugin.cpp
OBJECTS  := $(SOURCES:.cpp=.o)
TARGET   := aaOceanShaderDefs.so

$(TARGET) : $(OBJECTS)
	$(CXX) $(inputs) -shared $< -o $@ $(LDFLAGS)

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

all: $(TARGET)

clean:
	rm $(OBJECTS) $(TARGET)
