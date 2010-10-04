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
	$(OBJ)/aspicDataBasesPrinter.o 

############################################
############################################
LIBPATH = -L$(XERCESCROOT)/lib -L$(ASPICROOT)/lib

INCPATH = -I$(SRC) \
	-I$(XERCESCROOT)/include \
	-I$(ASPICROOT)/include

ifeq ($(MODE),Debug)
CXX_C_OPTIONS = -Wall -g
TARGET = $(BIN)/aspicDataBasesPrinterDebug.exe
LIBRARY = -lchemicsDebug -lxmlParserDebug -lgaussianDebug -lpolynomeDebug -lcontainorDebug -laspicUtilsDebug -lxerces-c 
endif

ifeq ($(MODE),Release)
CXX_C_OPTIONS = -O3 -DNDEBUG
TARGET = $(BIN)/aspicDataBasesPrinter.exe
LIBRARY = -lchemics -lxmlParser -lgaussian -lpolynome -lcontainor -laspicUtils -lxerces-c
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
