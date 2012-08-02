### Maria's Makefile ########################################################################################################################
###                                                                                                                                      ####
### gcc -o con my_conrec_exam_chris.c paulslib.c -lm bitmaplib.c -llapacke -llapack -lrefblas -L/usr/lib/gcc/x86_64-linux-gnu/4.4 -lgfortran ####
#############################################################################################################################################
OBJ  		= my_conrec_exam_chris.o 
LINKOBJ 	= paulslib.o bitmaplib.o 
BIN 		= my_conrec_exam_chris 

INC_PATH	= 
LIB_PATH 	= /usr/lib/gcc/x86_64-linux-gnu/4.4
LIBS 		= -lm -llapacke -llapack -lrefblas -lgfortran
CFLAGS 		= -I$(INC_PATH) -L$(LIB_PATH) 
RM = rm -f


### Actions ###
.PHONY: all

all: my_conrec_exam_chris

cleanObj:
	${RM} $(OBJ) $(LINKOBJ)

clean:
	${RM} $(OBJ) $(LINKOBJ) $(BIN)


### Executables ###
my_conrec_exam_chris: my_conrec_exam_chris.o $(LINKOBJ)
	$(CC) -o my_conrec_exam_chris my_conrec_exam_chris.o $(LINKOBJ) $(LIBS)


### Object files (for executables) ###
my_conrec_exam_chris.o: my_conrec_exam_chris.c
	$(CC) -c my_conrec_exam_chris.c -o my_conrec_exam_chris.o $(CFLAGS)


### Object files (for linkage) ###
paulslib.o: paulslib.c paulslib.h
	$(CC) -c  paulslib.c -o paulslib.o $(CFLAGS)

bitmaplib.o: bitmaplib.c bitmaplib.h
	$(CC) -c  bitmaplib.c -o bitmaplib.o $(CFLAGS)

