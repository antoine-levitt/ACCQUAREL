#########################################
#########################################
ifndef MODE 
MODE = Release
endif

ROOT = ../../..
INC = $(ROOT)/include
LIB = $(ROOT)/lib
OBJ = $(ROOT)/lib/obj/$(MODE)
SRC = ../src


OBJECTS = \
	$(OBJ)/aspicConfiguration.o \
	$(OBJ)/cmdLineParser.o
	
############################################
############################################
ifeq ($(MODE),Debug)
CXX_C_OPTIONS = -Wall -g
TARGET = $(LIB)/libaspicUtilsDebug.a
endif

ifeq ($(MODE),Release)
CXX_C_OPTIONS = -O3 -DNDEBUG
TARGET = $(LIB)/libaspicUtils.a
endif

INCPATH = -I$(SRC) -I../../libxmlParser/src -I../../libcontainor/src -I../../libpolynome/src -I../../libgaussian/src -I../../libchemics/src -I../../libgIntegrals/src -I$(XERCESCROOT)/include
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

exportIncludes :
	cp $(SRC)/*.h $(INC)

############################################
############################################
$(OBJ)/aspicConfiguration.o : $(SRC)/aspicConfiguration.cpp $(SRC)/aspicConfiguration.h

$(OBJ)/cmdLineParser.o : $(SRC)/cmdLineParser.cpp $(SRC)/cmdLineParser.h

$(OBJ)/mapleDebugInterface.o : $(SRC)/mapleDebugInterface.cpp $(SRC)/mapleDebugInterface.h

$(OBJ)/xmlBasisListParser.o : $(SRC)/xmlBasisListParser.cpp $(SRC)/xmlBasisListParser.h

