### Maria's Makefile ######################################
###                                                     ###
### Compiled commanmd:                                  ###
###                                                     ###
### gcc -o con                                          ###
###	conrec-example.c paulslib.c                     ###
###	-lm bitmaplib.c -llapacke -llapack -lrefblas    ###
###	-L/usr/lib/gcc/x86_64-linux-gnu/4.4 -lgfortran  ###
###########################################################

OBJ  		= conrec-example.o exclusion.o
LINKOBJ 	= paulslib.o bitmaplib.o 
BIN 		= conrec-example exclusion
INC_PATH	= 
LIB_PATH 	= /usr/lib/gcc/x86_64-linux-gnu/4.4
LIBS 		= -lm -llapacke -llapack -lrefblas -lgfortran
CFLAGS 		= -I$(INC_PATH) -L$(LIB_PATH) -Wall
RM = rm -f


### Actions ###
.PHONY: all

all: conrec-example exclusion

cleanObj:
	${RM} $(OBJ) $(LINKOBJ)

clean:
	${RM} $(OBJ) $(LINKOBJ) $(BIN)


### Executables ###
conrec-example: conrec-example.o $(LINKOBJ)
	$(CC) -o conrec-example conrec-example.o $(LINKOBJ) $(LIBS)

exclusion: exclusion.o $(LINKOBJ)
	$(CC) -o exclusion exclusion.o $(LINKOBJ) $(LIBS)


### Object files (src) ###
conrec-example.o: src/conrec-example.c
	$(CC) -c src/conrec-example.c -o conrec-example.o $(CFLAGS)

exclusion.o: src/exclusion.c
	$(CC) -c src/exclusion.c -o exclusion.o $(CFLAGS)


### Object files (lib) ###
paulslib.o: lib/paulslib.c lib/paulslib.h
	$(CC) -c  lib/paulslib.c -o paulslib.o $(CFLAGS)

bitmaplib.o: lib/bitmaplib.c lib/bitmaplib.h
	$(CC) -c  lib/bitmaplib.c -o bitmaplib.o $(CFLAGS)

