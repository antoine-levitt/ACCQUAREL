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
	$(OBJ)/centerPolynome.o \
	$(OBJ)/centerPolynome3D.o \
	$(OBJ)/complexPolynome.o \
	$(OBJ)/complexPolynome3D.o \
	$(OBJ)/polynome.o \
	$(OBJ)/polynome3D.o \
	$(OBJ)/polynomeBase.o \
	$(OBJ)/polynomeBase3D.o \
	$(OBJ)/string4Polynomes.o
	
############################################
############################################
ifeq ($(MODE),Debug)
CXX_C_OPTIONS = -Wall -g
TARGET = $(LIB)/libpolynomeDebug.a
endif

ifeq ($(MODE),Release)
CXX_C_OPTIONS = -O3 -DNDEBUG
TARGET = $(LIB)/libpolynome.a
endif

INCPATH = -I$(SRC) -I../../libcontainor/src -I../../libcomplex/src

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

$(OBJ)/centerPolynome.o : $(SRC)/centerPolynome.cpp $(SRC)/centerPolynome.h \
	$(SRC)/polynome.h

$(OBJ)/complexPolynome.o : $(SRC)/complexPolynome.cpp $(SRC)/complexPolynome.h \
	$(SRC)/polynomeBase.h

$(OBJ)/complexPolynome3D.o : $(SRC)/complexPolynome3D.cpp $(SRC)/complexPolynome3D.h \
	$(SRC)/polynomeBase3D.h

$(OBJ)/polynome.o : $(SRC)/polynome.cpp $(SRC)/polynome.h \
	$(SRC)/polynomeBase.h

$(OBJ)/polynome3D.o : $(SRC)/polynome3D.cpp $(SRC)/polynome3D.h \
	$(SRC)/polynomeBase3D.h

$(OBJ)/polynomeBase3D.o : $(SRC)/polynomeBase3D.cpp $(SRC)/polynomeBase3D.h \
	$(SRC)/string4Polynomes.h

$(OBJ)/polynomeBase.o : $(SRC)/polynomeBase.cpp $(SRC)/polynomeBase.h \
	$(SRC)/string4Polynomes.h

$(OBJ)/string4Polynomes.o : $(SRC)/string4Polynomes.cpp

