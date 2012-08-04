### Maria's Makefile ######################################
###                                                     ###
###########################################################

OBJ  		= exclusion_none.o exclusion_disk.o
LINKOBJ 	= paulslib.o bitmaplib.o 
BIN 		= exclusion_none exclusion_disk
INC_PATH	= 
LIB_PATH 	= /usr/lib/gcc/x86_64-linux-gnu/4.4
LIBS 		= -lm -llapacke -llapack -lrefblas -lgfortran
CFLAGS 		= -I$(INC_PATH) -L$(LIB_PATH) -Wall
RM = rm -f


### Actions ###
.PHONY: all

all: exclusion_none exclusion_disk

cleanObj:
	${RM} $(OBJ) $(LINKOBJ)

clean:
	${RM} $(OBJ) $(LINKOBJ) $(BIN)


### Executables ###
exclusion_none: exclusion_none.o $(LINKOBJ)
	$(CC) -o exclusion_none exclusion_none.o $(LINKOBJ) $(LIBS)

exclusion_disk: exclusion_disk.o $(LINKOBJ)
	$(CC) -o exclusion_disk exclusion_disk.o $(LINKOBJ) $(LIBS)


### Object files (src) ###
exclusion_none.o: src/exclusion_none.c
	$(CC) -c src/exclusion_none.c -o exclusion_none.o $(CFLAGS)

exclusion_disk.o: src/exclusion_disk.c
	$(CC) -c src/exclusion_disk.c -o exclusion_disk.o $(CFLAGS)


### Object files (lib) ###
paulslib.o: lib/paulslib.c lib/paulslib.h
	$(CC) -c  lib/paulslib.c -o paulslib.o $(CFLAGS)

bitmaplib.o: lib/bitmaplib.c lib/bitmaplib.h
	$(CC) -c  lib/bitmaplib.c -o bitmaplib.o $(CFLAGS)
	
