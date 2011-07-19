OBJ_DIR = ../../../obj/
BIN_DIR = ../../../bin/
# -------------------
# define our includes
# -------------------
INCLUDES = 

# ----------------------------------
# define our source and object files
# ----------------------------------
SOURCES=Align.cpp stdaln.c
OBJECTS=Align.o stdaln.o
_EXT_OBJECTS =
EXT_OBJECTS=$(patsubst %,$(OBJ_DIR)/%,$(_EXT_OBJECTS))
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

$(BUILT_OBJECTS): $(SOURCES)
	gcc -Wall -Werror -O3 -c stdaln.c  -o $(OBJ_DIR)/stdaln.o
	g++ -Wall -Werror -O3 -c Align.cpp -o $(OBJ_DIR)/Align.o
	

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/* $(BIN_DIR)/*

.PHONY: clean