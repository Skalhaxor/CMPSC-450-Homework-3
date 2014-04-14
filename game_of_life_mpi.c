// Created by: Alexander Anderson
// CMPSC 450 (Concurrent Programming)
// Game of Life (MPI)
// Date: 4/11/2014

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "mpi.h"

#define GRID_SIZE  30                // n by n size of grid
#define TIME_STEPS 100               // number of times to update
#define SLEEP_TIME 250 * 1000        // millisec * microsec-conversion of time to sleep

#define TOP          0
#define BOTTOM       1
#define TOP_RIGHT    2
#define TOP_LEFT     3
#define BOTTOM_RIGHT 4
#define BOTTOM_LEFT  5
#define COUNT        6

int** AllocateGrid(int);
void  DeallocateGrid(int**, int);
void  InitGrid(int**, int);
int   CountLive(int**, int);

int main(int argc, char* argv[])
{
    int** grid;             // grid containing the cells and whether they are alive or dead
    int myRank;             // id for each processor
    int aboveProc;          // id of process above the current (proc 0 and last one loop around)
    int belowProc;          // id of process below the current (proc 0 and last one loop around)
    int numProcs;           // number of processors
    int numRows;            // row size for the grid for each processor
    int numLive;            // number of cells alive
    MPI_Status status;

    // initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    if(myRank == 0)
        aboveProc = numProcs - 1;
    else
        aboveProc = myRank - 1;
    
    belowProc = (myRank+1) % numProcs;

    // calculate the number of rows the grid for each processor will have
    if(myRank == numProcs - 1)      // last processor
        numRows = (GRID_SIZE / numProcs) + (GRID_SIZE % numProcs);
    else                            // all other processors
        numRows = GRID_SIZE / numProcs;
        
    // allocate the grid
    grid = AllocateGrid(numRows);
    if(grid == NULL)
    {
        printf("ERROR: Could not allocate memory for grid\n");
        return -1;
    }
    
    // initialize and print the grid
    InitGrid(grid, numRows);
    
    // perform the game of life for TIME_STEPS iterations
    for(int time = 1; time <= TIME_STEPS; ++time)
    {
        // sync left and right boundaries of grid
        for(int i = 1; i <= numRows; ++i)
        {
            grid[i][0] = grid[i][GRID_SIZE];
            grid[i][GRID_SIZE+1] = grid[i][1];
        }
        
        // send top boundary of grid
        MPI_Sendrecv(&grid[1][1],         GRID_SIZE, MPI_INT, aboveProc, TOP,
                     &grid[numRows+1][1], GRID_SIZE, MPI_INT, belowProc, TOP,
                                                      MPI_COMM_WORLD, &status);
        
        // send bottom boundary of grid
        MPI_Sendrecv(&grid[numRows][1], GRID_SIZE, MPI_INT, belowProc, BOTTOM,
                     &grid[0][1],       GRID_SIZE, MPI_INT, aboveProc, BOTTOM,
                                                       MPI_COMM_WORLD, &status);
        
        // send top right corner of grid
        MPI_Sendrecv(&grid[1][GRID_SIZE], 1, MPI_INT, aboveProc, TOP_RIGHT,
                     &grid[numRows+1][0], 1, MPI_INT, belowProc, TOP_RIGHT,
                                                    MPI_COMM_WORLD, &status);
        
        // send top left corner of grid
        MPI_Sendrecv(&grid[1][1],                   1, MPI_INT, aboveProc, TOP_LEFT,
                     &grid[numRows+1][GRID_SIZE+1], 1, MPI_INT, belowProc, TOP_LEFT,
                                                             MPI_COMM_WORLD, &status);
        
        // send bottom right corner of grid
        MPI_Sendrecv(&grid[numRows][GRID_SIZE], 1, MPI_INT, belowProc, BOTTOM_RIGHT,
                     &grid[0][0],               1, MPI_INT, aboveProc, BOTTOM_RIGHT,
                                                             MPI_COMM_WORLD, &status);
        
        // send bottom left corner of grid
        MPI_Sendrecv(&grid[numRows][1],     1, MPI_INT, belowProc, BOTTOM_LEFT,
                     &grid[0][GRID_SIZE+1], 1, MPI_INT, aboveProc, BOTTOM_LEFT,
                                                        MPI_COMM_WORLD, &status);
    
        // allocate temporary grid for storing the simultaneous update
        int** tempGrid = AllocateGrid(numRows);
        if(tempGrid == NULL)
        {
            printf("ERROR: Could not allocate memory for grid\n");
            return -1;
        }
    
        for(int i = 1; i <= numRows; ++i)
        {
            for(int j = 1; j <= GRID_SIZE; ++j)
            {
                int numNeighborAlive = 0;
            
                for(int k = -1; k <= 1; ++k)
                {
                    for(int l = -1; l <= 1; ++l)
                    {
                        int row = i + k;
                        int col = j + l;
                    
                        if( (row != i) || (col != j) )
                            numNeighborAlive += grid[row][col];
                    }
                }
                
                switch(numNeighborAlive)
                {
                case 2:
                    tempGrid[i][j] = grid[i][j];
                    break;
                    
                case 3:
                    tempGrid[i][j] = 1;
                    break;
                    
                default:
                    tempGrid[i][j] = 0;
                }
            }
        }
        
        DeallocateGrid(grid, numRows);
        grid = tempGrid;
    }
    
    numLive = CountLive(grid, numRows);
    
    if(myRank != 0)
    {
        MPI_Send(&numLive, 1, MPI_INT, 0, COUNT, MPI_COMM_WORLD);
    }
    else
    {
        int temp;
    
        for(int i = 1; i < numProcs; ++i)
        {
            MPI_Recv(&temp, 1, MPI_INT, i, COUNT, MPI_COMM_WORLD, &status);
            numLive += temp;
        }
        
        printf("Number of live cells = %d\n", numLive);
    }
    
    MPI_Finalize();
    DeallocateGrid(grid, numRows);

    return 0;
}

int** AllocateGrid(int numRows)
{
    int** grid;
    
    grid = malloc(sizeof(int) * (numRows+2));
    
    if(grid == NULL)
    {
        printf("ERROR: Could not allocate memory for grid\n");
        return NULL;
    }
    
    for(int i = 0; i < numRows+2; ++i)
    {
        grid[i] = malloc(sizeof(int) * (GRID_SIZE+2));
        
        if(grid == NULL)
        {
            printf("ERROR: Could not allocate memory for grid\n");
            return NULL;
        }
    }
    
    return grid;
}

void DeallocateGrid(int** grid, int numRows)
{
    for(int i = 0; i < numRows+2; ++i)
    {
        free(grid[i]);
    }

    free(grid);
}

void InitGrid(int** grid, int numRows)
{
    for(int i = 1; i <= numRows; ++i)
    {
        for(int j = 1; j <= GRID_SIZE; ++j)
        {
            grid[i][j] = rand() % 2;
        }
    }
}

int CountLive(int** grid, int numRows)
{
    int numLive = 0;

    for(int i = 1; i <= numRows; ++i)
    {
        for(int j = 1; j <= GRID_SIZE; ++j)
        {
            numLive += grid[i][j];
        }
    }
    
    return numLive;
}

