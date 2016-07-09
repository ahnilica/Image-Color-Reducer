CC = g++
AR = ar
INCLUDES = -I../include
CFLAGS = $(INCLUDES)
CPP_FILES := $(wildcard *.cpp)
OBJ_FILES := $(notdir $(CPP_FILES:.cpp=.o))

.PHONY: all clean

all: libmatrix.a

libmatrix.a: $(OBJ_FILES)
	$(AR) rvu $@ $(OBJ_FILES)
	ranlib $@

clean:
	rm -f *.o libmatrix.a

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@
