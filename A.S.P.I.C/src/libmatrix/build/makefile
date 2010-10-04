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
	$(OBJ)/generalizedEigenSolver.o \
	$(OBJ)/matrix.o \
	$(OBJ)/matrixFull.o \
	$(OBJ)/matrixNonSparse.o \
	$(OBJ)/matrixSymetric.o \
	$(OBJ)/tensor.o \
	$(OBJ)/tensorFull.o \
	$(OBJ)/tensorNonSparse.o \
	$(OBJ)/tensorSymetric.o \
	$(OBJ)/vector.o \
	$(OBJ)/xmlMatrixDocumentParser.o \
	$(OBJ)/xmlMatrixParser.o \
	$(OBJ)/xmlTensorDocumentParser.o \
	$(OBJ)/xmlTensorParser.o
	
############################################
############################################
ifeq ($(MODE),Debug)
CXX_C_OPTIONS = -Wall -g
TARGET = $(LIB)/libmatrixDebug.a
endif

ifeq ($(MODE),Release)
CXX_C_OPTIONS = -O3 -DNDEBUG
TARGET = $(LIB)/libmatrix.a
endif

INCPATH = -I$(SRC)  \
	-I../../libcontainor/src \
	-I../../libxmlParser/src \
	-I../../aspicUtils/src \
	-I$(XERCESCROOT)/include

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
