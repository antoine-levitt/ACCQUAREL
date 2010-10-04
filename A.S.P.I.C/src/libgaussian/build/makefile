#########################################
#########################################
ifndef MODE 
MODE = Release
endif

OS = macOsX

ROOT = ../../..
INC = $(ROOT)/include
LIB = $(ROOT)/lib
OBJ = $(ROOT)/lib/obj/$(MODE)
SRC = ../src

OBJECTS = \
	$(OBJ)/error_function.o \
	$(OBJ)/gaussian.o \
	$(OBJ)/gaussian3D.o \
	$(OBJ)/gaussianPolynome.o \
	$(OBJ)/gaussianPolynome3D.o 
	
############################################
############################################
ifeq ($(MODE),Debug)
CXX_C_OPTIONS = -Wall -g
TARGET = $(LIB)/libgaussianDebug.a
endif

ifeq ($(MODE),Release)
CXX_C_OPTIONS = -O3 -DNDEBUG
TARGET = $(LIB)/libgaussian.a
endif

INCPATH = -I$(SRC) -I../../libcontainor/src -I../../libpolynome/src

LIBPATH = 

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
$(OBJ)/gaussian.o : $(SRC)/gaussian.cpp $(SRC)/gaussian.h

$(OBJ)/gaussian_3D.o : $(SRC)/gaussian_3D.cpp $(SRC)/gaussian_3D.h

$(OBJ)/gaussian_polynome.o : $(SRC)/gaussian_polynome.cpp $(SRC)/gaussian_polynome.h \
	$(SRC)/gaussian.h

$(OBJ)/gaussian_polynome_3D.o : $(SRC)/gaussian_polynome.cpp $(SRC)/gaussian_polynome.h \
	$(SRC)/gaussian_3D.h


