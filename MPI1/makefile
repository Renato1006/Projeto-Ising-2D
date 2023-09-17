IDIR=.
CC=mpicc
CFLAGS=-I$(IDIR) -O3 -fexpensive-optimizations -m64 -foptimize-register-move -funroll-loops -ffast-math -mtune=native -march=native

ODIR=obj
LDIR=../lib
BIN=bin
SRC=src

_DEPS=ising_paralel1.h
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ=ising_paralel1.o
OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ))

LIBS=-lm

$(ODIR)/%.o: $(SRC)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

ising.x: $(OBJ)
	$(CC) -o $(BIN)/$@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -f $(ODIR)/*.o $(BIN)/*.x *~ core 
