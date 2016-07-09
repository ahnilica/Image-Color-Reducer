CC = g++
INCLUDES = -I./Matrix/include
LIBS = -L./Matrix/lib -lmatrix
CFLAGS = $(INCLUDES)

.PHONY: all clean

all: WTA

WTA: WTA.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean :
	-rm *.o

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@
