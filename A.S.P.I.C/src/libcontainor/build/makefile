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
	$(OBJ)/containor.o \
	$(OBJ)/containor2D.o \
	$(OBJ)/containor3D.o \
	$(OBJ)/containorSparse.o \
	$(OBJ)/containorSparse2D.o \
	$(OBJ)/containorSparse3D.o \
	$(OBJ)/containorTri2D.o \
	$(OBJ)/dpoint.o \
	$(OBJ)/ipoint.o \
	$(OBJ)/tpoint.o \
	$(OBJ)/xmlDPoint3DocumentParser.o \
	$(OBJ)/xmlDPoint3Parser.o \
	$(OBJ)/xmlIPoint3Parser.o \
	$(OBJ)/xmlIPoint3DocumentParser.o

############################################
############################################
ifeq ($(MODE),Debug)
TARGET = $(LIB)/libcontainorDebug.a
CXX_C_OPTIONS = -Wall -g
endif

ifeq ($(MODE),Release)
TARGET = $(LIB)/libcontainor.a
CXX_C_OPTIONS = -O3 -DNDEBUG
endif

INCPATH = -I$(SRC) \
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
$(OBJ)/containor.o : $(SRC)/containor.cpp $(SRC)/containor.h

$(OBJ)/containor2D.o : $(SRC)/containor2D.cpp $(SRC)/containor2D.h \
	$(SRC)/containor.h

$(OBJ)/containor3D.o : $(SRC)/containor3D.cpp $(SRC)/containor3D.h \
	$(SRC)/containor.h

$(OBJ)/containorSparse.o : $(SRC)/containorSparse.cpp $(SRC)/containorSparse.h

$(OBJ)/containorSparse2D.o : $(SRC)/containorSparse2D.cpp $(SRC)/containorSparse2D.h \
	$(SRC)/ipoint.h

$(OBJ)/containorSparse3D.o : $(SRC)/containorSparse3D.cpp $(SRC)/containorSparse3D.h \
	$(SRC)/ipoint.h

$(OBJ)/containorTri2D.o : $(SRC)/containorTri2D.cpp $(SRC)/containorTri2D.h \
	$(SRC)/containor.h

$(OBJ)/tpoint.o : $(SRC)/tpoint.cpp $(SRC)/tpoint.h

$(OBJ)/dpoint.o : $(SRC)/dpoint.cpp $(SRC)/dpoint.h \
	$(SRC)/tpoint.h

$(OBJ)/ipoint.o : $(SRC)/ipoint.cpp $(SRC)/ipoint.h \
	$(SRC)/tpoint.h

$(OBJ)/xmlDPoint3DocumentParser.o :	 $(SRC)/xmlDPoint3DocumentParser.cpp $(SRC)/xmlDPoint3DocumentParser.h \
	$(SRC)/xmlDPoint3Parser.h

$(OBJ)/xmlDPoint3Parser.o : $(SRC)/xmlDPoint3Parser.cpp $(SRC)/xmlDPoint3Parser.h


$(OBJ)/xmlIPoint3Parser.o : $(SRC)/xmlIPoint3Parser.cpp $(SRC)/xmlIPoint3Parser.h

############################################
