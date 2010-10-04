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
	$(OBJ)/xmlDataBaseInterface.o \
	$(OBJ)/xmlDocumentParser.o \
	$(OBJ)/xmlErrorReporter.o \
	$(OBJ)/xmlParser.o

############################################
############################################
ifeq ($(MODE),Debug)
TARGET = $(LIB)/libxmlParserDebug.a
CXX_C_OPTIONS = -Wall -g
endif

ifeq ($(MODE),Release)
TARGET = $(LIB)/libxmlParser.a
CXX_C_OPTIONS = -O3 -DNDEBUG
endif

INCPATH = -I$(SRC) -I$(XERCESCROOT)/include

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
$(OBJ)/xmlDataBaseInterface.o : $(SRC)/xmlDataBaseInterface.cpp $(SRC)/xmlDataBaseInterface.h \
	$(SRC)/xmlDocumentParser.h

$(OBJ)/xmlDocumentParser.o : $(SRC)/xmlDocumentParser.cpp $(SRC)/xmlDocumentParser.h \
	$(SRC)/xmlErrorReporter.h \
	$(SRC)/xmlParser.h

$(OBJ)/xmlErrorReporter.o : $(SRC)/xmlErrorReporter.cpp $(SRC)/xmlErrorReporter.h \
	$(SRC)/xmlErrorReporter.h

$(OBJ)/xmlParser.o : $(SRC)/xmlParser.cpp $(SRC)/xmlParser.h \
	$(SRC)/xmlErrorReporter.h
