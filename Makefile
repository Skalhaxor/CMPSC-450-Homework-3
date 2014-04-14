CC=gcc
MPICC=mpicc
CFLAGS=-std=c99 -Wall -Wextra
LDFLAGS=-lncurses
SERIAL_NAME=gol-serial
SERIAL_FILE=game_of_life_serial.c
MPI_NAME=gol-mpi
MPI_FILE=game_of_life_mpi.c

serial:
	$(CC) $(CFLAGS) -o $(SERIAL_NAME) $(SERIAL_FILE) $(LDFLAGS)

serial-dbg:
	$(CC) $(CFLAGS) -g -o $(SERIAL_NAME) $(SERIAL_FILE) $(LDFLAGS)

mpi:
	$(MPICC) $(CFLAGS) -o $(MPI_NAME) $(MPI_FILE) $(LDFLAGS)

mpi-dbg:
	$(MPICC) $(CFLAGS) -g -o $(MPI_NAME) $(MPI_FILE) $(LDFLAGS)

