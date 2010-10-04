#########################################
#########################################
ifndef MODE 
MODE = Release
endif

ROOT = ../../..

INC = $(ROOT)/include
LIB = $(ROOT)/lib
BIN = $(ROOT)/bin
OBJ = $(ROOT)/lib/obj/$(MODE)
SRC = ../src

OBJECTS = \
	$(OBJ)/coulombIntegral.o \
	$(OBJ)/hermiteIntegrationBuffer.o \
	$(OBJ)/hermiteProjector.o \
	$(OBJ)/monoElectronicIntegral.o \
	$(OBJ)/potentialIntegral.o

############################################
############################################
ifeq ($(MODE),Debug)
TARGET = $(LIB)/libgIntegralsDebug.a
CXX_C_OPTIONS = -Wall -g
endif

ifeq ($(MODE),Release)
TARGET = $(LIB)/libgIntegrals.a
CXX_C_OPTIONS = -O3 -DNDEBUG
endif

INCPATH = -I$(SRC) -I../../libxmlParser/src -I$(XERCESCROOT)/include -I../../libchemics/src -I../../libpolynome/src -I../../libgaussian/src -I../../libcontainor/src

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
