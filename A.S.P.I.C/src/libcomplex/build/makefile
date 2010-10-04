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
	$(OBJ)/complex.o 
	
############################################
############################################
ifeq ($(MODE),Debug)
CXX_C_OPTIONS = -Wall -g
TARGET = $(LIB)/libcomplexDebug.a
endif

ifeq ($(MODE),Release)
CXX_C_OPTIONS = -O3 -DNDEBUG
TARGET = $(LIB)/libcomplex.a
endif

INCPATH = -I$(SRC) 
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
$(OBJ)/complex.o : $(SRC)/complex.cpp $(SRC)/complex.h
