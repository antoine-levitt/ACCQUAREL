#########################################
#########################################
ifndef MODE 
MODE = Release
endif

#########################################
#########################################
ROOT = ../../..
INC = $(ROOT)/include
LIB = $(ROOT)/lib
OBJ = $(ROOT)/lib/obj/$(MODE)
BIN = $(ROOT)/bin
SRC = ../src

OBJECTS_OVERLAP = \
	$(OBJ)/overlapIntegralSample.o \
	$(OBJ)/mapleDebugInterface.o \

OBJECTS_KINETIC = \
	$(OBJ)/kineticIntegralSample.o \
	$(OBJ)/mapleDebugInterface.o \

OBJECTS_COULOMB = \
	$(OBJ)/coulombIntegralSample.o \
	$(OBJ)/mapleDebugInterface.o \

OBJECTS_POTENTIAL = \
	$(OBJ)/potentialIntegralSample.o \
	$(OBJ)/mapleDebugInterface.o \

############################################
############################################
LIBPATH = -L$(XERCESCROOT)/lib -L$(ASPICROOT)/lib

INCPATH = -I$(SRC) \
	-I$(XERCESCROOT)/include \
	-I$(ASPICROOT)/include

ifeq ($(MODE),Debug)
CXX_C_OPTIONS = -Wall -g
TARGET_OVERLAP = $(BIN)/overlapIntegralSampleDebug.exe
TARGET_KINETIC = $(BIN)/kineticIntegralSampleDebug.exe
TARGET_COULOMB = $(BIN)/coulombIntegralSampleDebug.exe
TARGET_POTENTIAL= $(BIN)/potentialIntegralSampleDebug.exe
LIBRARY = -lgIntegralsDebug -lchemicsDebug -lxmlParserDebug -lgaussianDebug -lpolynomeDebug -lcontainorDebug -laspicUtilsDebug -lxerces-c 
endif

ifeq ($(MODE),Release)
CXX_C_OPTIONS = -O3 -DNDEBUG
TARGET_OVERLAP = $(BIN)/overlapIntegralSample.exe
TARGET_KINETIC = $(BIN)/kineticIntegralSample.exe
TARGET_COULOMB = $(BIN)/coulombIntegralSample.exe
TARGET_POTENTIAL = $(BIN)/potentialIntegralSample.exe
LIBRARY = -lgIntegrals -lchemics -lxmlParser -lgaussian -lpolynome -lcontainor -laspicUtils -lxerces-c
endif

############################################
############################################
$(LIB)/%.a :
	ar -rcuv $@ $^

$(BIN)/%.exe :
	$(CXX) -o $@ $^ $(LIBPATH) $(LIBRARY)

$(OBJ)/%.o : $(SRC)/%.cpp
	$(CXX) -c $(CXX_C_OPTIONS) -o $@ $(INCPATH) $<

############################################
############################################
all : $(TARGET_OVERLAP) $(TARGET_KINETIC) $(TARGET_COULOMB) $(TARGET_POTENTIAL)

clean :
	-rm -f $(OBJECTS_OVERLAP) $(OBJECTS_KINETIC) $(OBJECTS_COULOMB) $(OBJECTS_POTENTIAL)
	-rm -f $(TARGET_OVERLAP) $(TARGET_KINETIC) $(TARGET_COULOMB) $(TARGET_POTENTIAL) 
	
$(TARGET_OVERLAP) : $(OBJECTS_OVERLAP)

$(TARGET_KINETIC) : $(OBJECTS_KINETIC)

$(TARGET_COULOMB) : $(OBJECTS_COULOMB)

$(TARGET_POTENTIAL) : $(OBJECTS_POTENTIAL)
#############################################
