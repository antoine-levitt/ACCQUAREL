ifndef MODE 
MODE = Release
endif

all : exportIncludes
	cd aspicUtils/build && MODE=$(MODE) make
	cd libchemics/build && MODE=$(MODE) make
	cd libcomplex/build && MODE=$(MODE) make
	cd libcontainor/build && MODE=$(MODE) make
	cd libgaussian/build && MODE=$(MODE) make
	cd libgIntegrals/build && MODE=$(MODE) make
	cd libmatrix/build && MODE=$(MODE) make
	cd libmodel/build && MODE=$(MODE) make
	cd liboptimisation/build && MODE=$(MODE) make
	cd libpolynome/build && MODE=$(MODE) make
	cd libxmlParser/build && MODE=$(MODE) make

clean :
	cd aspicUtils/build && MODE=$(MODE) make clean
	cd libchemics/build && MODE=$(MODE) make clean
	cd libcomplex/build && MODE=$(MODE) make clean
	cd libcontainor/build && MODE=$(MODE) make clean
	cd libgaussian/build && MODE=$(MODE) make clean
	cd libgIntegrals/build && MODE=$(MODE) make clean
	cd libmatrix/build && MODE=$(MODE) make clean
	cd libmodel/build && MODE=$(MODE) make clean
	cd liboptimisation/build && MODE=$(MODE) make clean
	cd libpolynome/build && MODE=$(MODE) make clean
	cd libxmlParser/build && MODE=$(MODE) make clean

exportIncludes :
	cd aspicUtils/build && MODE=$(MODE) make exportIncludes
	cd libchemics/build && MODE=$(MODE) make exportIncludes
	cd libcomplex/build && MODE=$(MODE) make exportIncludes
	cd libcontainor/build && MODE=$(MODE) make exportIncludes
	cd libgaussian/build && MODE=$(MODE) make exportIncludes
	cd libgIntegrals/build && MODE=$(MODE) make exportIncludes
	cd libmatrix/build && MODE=$(MODE) make exportIncludes
	cd libmodel/build && MODE=$(MODE) make exportIncludes
	cd liboptimisation/build && MODE=$(MODE) make exportIncludes
	cd libpolynome/build && MODE=$(MODE) make exportIncludes
	cd libxmlParser/build && MODE=$(MODE) make exportIncludes
