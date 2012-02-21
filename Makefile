CC=gcc
FLAGS=-Wall -lm
EXEC=wave
EXEC_TEST=wave_test
CFILES=main.c wave.c

all: ${CFILES}
	${CC} ${FLAGS} -o ${EXEC} ${CFILES} -fopenmp


pwave: print_wave.o wave.o spline.o
	${CC} ${FLAGS} -o pwave print_wave.o wave.o 

seq:
	${CC} ${FLAGS} -o ${EXEC}_seq ${CFILES}

clean:
	rm -f *.o *.swp ${EXEC}
