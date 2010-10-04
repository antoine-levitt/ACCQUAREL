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

OBJECTS = \
	$(OBJ)/generalizedEigenSolverSample.o 

############################################
############################################
LIBPATH = -L$(XERCESCROOT)/lib -L$(ASPICROOT)/lib -L/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A


INCPATH = -I$(SRC) \
	-I$(XERCESCROOT)/include \
	-I$(ASPICROOT)/include

ifeq ($(MODE),Debug)
CXX_C_OPTIONS = -Wall -g
TARGET = $(BIN)/generalizedEigenSolverSampleDebug.exe
LIBRARY = -lmatrixDebug -lxmlParserDebug -lcontainorDebug -laspicUtilsDebug -lxerces-c -llapack -lblas
endif

ifeq ($(MODE),Release)
CXX_C_OPTIONS = -O3 -DNDEBUG
TARGET = $(BIN)/generalizedEigenSolverSample.exe
LIBRARY = -lmatrix -lxmlParser -lcontainor -laspicUtils -lxerces-c  -llapack -lblas 
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
all : $(TARGET)

clean :
	-rm -f $(OBJECTS)
	-rm -f $(TARGET)
	
$(TARGET) : $(OBJECTS)
