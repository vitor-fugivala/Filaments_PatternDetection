CC=g++ -fopenmp
INC_PATH = ./Eigen
CFLAGS= -c -Wall -Werror -I$(INC_PATH)
LDFLAGS=
SOURCES= kdtree_2d.hpp kdtree_2d.cpp find_fil_final.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=find_filament

all: $(SOURCES) $(EXECUTABLE)
    
$(EXECUTABLE): $(OBJECTS) ;$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o: ;$(CC) $(CFLAGS) $< -o $@


