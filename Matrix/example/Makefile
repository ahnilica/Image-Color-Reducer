CC = g++
INCLUDES = -I../include
LIBS = -L../lib -lmatrix
CFLAGS = $(INCLUDES)

EXECUTABLES = testfileIO testimageIO testMatrix testEig testMpp proj1

.PHONY: all clean

all: $(EXECUTABLES)

proj1: proj1.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

testMpp: testMpp.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

testEig: testEig.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

testfileIO: testfileIO.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

testimageIO: testimageIO.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

testMatrix: testMatrix.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(EXECUTABLES)

%.o: %.cpp
	$(CC) -c $(CFLAGS) $< -o $@
