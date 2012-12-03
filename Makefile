
INCLUDE_DIR = ../ptex/install/include
LIB_DIR = ../ptex/install/lib

LIBS = -lPtex

CFLAGS = $(DEFINES) -c -g -I$(INCLUDE_DIR)
LDFLAGS = -L$(LIB_DIR) -Wl,-rpath=$(LIB_DIR)

SOURCES = main.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all : $(SOURCES) ptxmerge

ptxmerge : $(OBJECTS)
	$(CXX) $(LDFLAGS) $(OBJECTS) -o $@ $(LIBS)

.cpp.o :
	$(CXX) $(CFLAGS) $< -o $@
