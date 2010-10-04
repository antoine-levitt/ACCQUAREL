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
SRC = ../src

OBJECTS = \
	$(OBJ)/atom.o \
	$(OBJ)/basisElement.o \
	$(OBJ)/chemicalElement.o \
	$(OBJ)/chemicalElementsDataBaseInterface.o \
	$(OBJ)/contractions.o \
	$(OBJ)/gaussianBasisFunction.o \
	$(OBJ)/gaussianBasisShell.o \
	$(OBJ)/molecule.o \
	$(OBJ)/moleculeMapper.o \
	$(OBJ)/shellType.o \
	$(OBJ)/xmlAtomParser.o \
	$(OBJ)/xmlBasisDataBaseInterface.o \
	$(OBJ)/xmlBasisElementParser.o \
	$(OBJ)/xmlBasisFunctionDocumentParser.o \
	$(OBJ)/xmlBasisFunctionParser.o \
	$(OBJ)/xmlBasisKeyListParser.o \
	$(OBJ)/xmlContractionsParser.o \
	$(OBJ)/xmlMoleculeDocumentParser.o \
	$(OBJ)/xmlMoleculeParser.o \
	$(OBJ)/xmlShellTypeParser.o \
	$(OBJ)/xmlShellTypesDataBaseInterface.o

############################################
############################################
ifeq ($(MODE),Debug)
CXX_C_OPTIONS = -Wall -g
TARGET = $(LIB)/libchemicsDebug.a
endif

ifeq ($(MODE),Release)
CXX_C_OPTIONS = -O3 -DNDEBUG
TARGET = $(LIB)/libchemics.a
endif

INCPATH = -I$(SRC) \
	-I$(XERCESCROOT)/include \
	-I../../aspicUtils/src \
	-I../../libxmlParser/src \
	-I../../libgaussian/src \
	-I../../libpolynome/src \
	-I../../libcontainor/src
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

## La classe atome hérite de la classe chemicalElement.
$(OBJ)/atom.o : $(SRC)/atom.cpp $(SRC)/atom.h \
	$(SRC)/chemicalElement.h

## La classe gaussianBasisElement
## dépend de gaussianBasisShell
## dépend de atom et donc de chemicalElement.
$(OBJ)/basisElement.o : $(SRC)/basisElement.cpp $(SRC)/basisElement.h \
	$(SRC)/contractions.h \
	$(SRC)/shellType.h 

## La classe chemical Element est une classe de base.
$(OBJ)/chemicalElement : $(SRC)/chemicalElement.cpp $(SRC)/chemicalElement.h

## la classe chemicalElementsDataBaseInterface :
## elle dépend de chemicalElement et de xmlChemicalElementParser.
$(OBJ)/chemicalElementsDataBaseInterface.o : $(SRC)/chemicalElementsDataBaseInterface.cpp $(SRC)/chemicalElementsDataBaseInterface.h \
	$(SRC)/chemicalElement.h

## La classe contraction ne dépend de rien.
$(OBJ)/contractions.o : $(SRC)/contractions.cpp $(SRC)/contractions.h

## La classe gaussianBasisFunction ne dépend que de contraction.
$(OBJ)/gaussianBasisFunction.o : $(SRC)/gaussianBasisFunction.cpp $(SRC)/gaussianBasisFunction.h \
	$(SRC)/contractions.h

## La classe gaussianBasisShell ne dépend que de contraction.
$(OBJ)/gaussianBasisShell.o : $(SRC)/gaussianBasisShell.cpp $(SRC)/gaussianBasisShell.h \
	$(SRC)/contractions.h \
	$(SRC)/gaussianBasisFunction.h \
	$(SRC)/shellType.h \

## la classe molécule dépend de atom
## donc de tout ce qui dépend d'atome.
$(OBJ)/molecule.o : $(SRC)/molecule.cpp $(SRC)/molecule.h \
	$(SRC)/atom.h \
	$(SRC)/chemicalElement.h \
	$(SRC)/basisElement.h \
	$(SRC)/gaussianBasisFunction.h \
	$(SRC)/gaussianBasisShell.h

## la classe gaussianBasisFunctionParser
## dépend de gaussianBasisFunction.hpp
## de contraction.hpp
$(OBJ)/moleculeMapper.o : $(SRC)/moleculeMapper.cpp $(SRC)/moleculeMapper.h \
	$(SRC)/contractions.h \
	$(SRC)/gaussianBasisFunction.h \
	$(SRC)/gaussianBasisShell.h \
	$(SRC)/shellType.h \
	$(SRC)/molecule.h

##
$(OBJ)/shellType.o : $(SRC)/shellType.cpp $(SRC)/shellType.h \
	$(SRC)/xmlShellTypesDataBaseInterface.h

## la classe 
$(OBJ)/xmlAtomParser.o : $(SRC)/xmlAtomParser.cpp $(SRC)/xmlAtomParser.h \
	$(SRC)/atom.h \
	$(SRC)/chemicalElement.h \
	$(SRC)/chemicalElementsDataBaseInterface.h

## La classe gaussianBasisDataBaseInterface :
## gaussianBasisElement
## xmlGaussianBasisElementParser
$(OBJ)/xmlBasisDataBaseInterface.o : $(SRC)/xmlBasisDataBaseInterface.cpp $(SRC)/xmlBasisDataBaseInterface.h \
	$(SRC)/xmlBasisElementParser.h

##
$(OBJ)/xmlContractionsParser.o : $(SRC)/xmlContractionsParser.cpp $(SRC)/xmlContractionsParser.h \
	$(SRC)/contractions.h

## la classe gaussianBasisFunctionParser
## dépend de gaussianBasisFunction.hpp
## de contraction.hpp
$(OBJ)/xmlBasisFunctionDocumentParser.o : $(SRC)/xmlBasisFunctionDocumentParser.cpp \
	$(SRC)/xmlBasisFunctionDocumentParser.h \
	$(SRC)/contractions.h \
	$(SRC)/gaussianBasisFunction.h \
	$(SRC)/xmlContractionsParser.h \
	$(SRC)/xmlBasisFunctionParser.h

## la classe gaussianBasisFunctionParser
## dépend de gaussianBasisFunction.hpp
## de contraction.hpp
$(OBJ)/xmlBasisFunctionParser.o : $(SRC)/xmlBasisFunctionParser.cpp \
	$(SRC)/xmlBasisFunctionParser.h \
	$(SRC)/contractions.h \
	$(SRC)/gaussianBasisFunction.h \
	$(SRC)/xmlContractionsParser.h

##
$(OBJ)/xmlMoleculeParser.o : $(SRC)/xmlMoleculeParser.cpp $(SRC)/xmlMoleculeParser.h \
	$(SRC)/atom.h \
	$(SRC)/xmlBasisDataBaseInterface.h \
	$(SRC)/basisElement.h \
	$(SRC)/molecule.h \
	$(SRC)/xmlAtomParser.h

##
$(OBJ)/xmlMoleculeDocumentParser.o : $(SRC)/xmlMoleculeDocumentParser.cpp $(SRC)/xmlMoleculeDocumentParser.h \
	$(SRC)/xmlMoleculeParser.h

##
$(SRC)/xmlShellTypeParser.o : $(SRC)/xmlShellTypeParser.cpp $(SRC)/xmlShellTypeParser.h \
	$(SRC)/shellType.h

##
$(OBJ)/xmlShellTypesDataBaseInterface.o : $(SRC)/xmlShellTypesDataBaseInterface.cpp $(SRC)/xmlShellTypesDataBaseInterface.h \
	$(SRC)/shellType.h \
	$(SRC)/xmlShellTypeParser.h 

