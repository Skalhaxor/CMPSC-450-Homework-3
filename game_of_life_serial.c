// Created by: Alexander Anderson
// CMPSC 450 (Concurrent Programming)
// Game of Life (Serial)
// Date: 4/9/2014

// define if you want text/graphical display of the grid that will update (be replaced) each iteration
//#define USE_CURSES

// define if you want the grid to be displayed
#define PRINT_GRID

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#ifdef USE_CURSES
#include <curses.h>
#endif

#define GRID_SIZE  30                // n by n size of grid
#define TIME_STEPS 100               // number of times to update
#define SLEEP_TIME 250 * 1000        // millisec * microsec-conversion of time to sleep

int** AllocateGrid();
void InitGrid(int**);
void PrintGrid(int**);
void CountLive(int**);

int main(int argc, char* argv[])
{
    int** grid;            // grid containing the cells and whether they are alive or dead
    
    // set up curses for displaying the grid
#ifdef USE_CURSES
    initscr();
    if(has_colors())
    {
        start_color();
        init_pair(1, COLOR_BLACK, COLOR_WHITE);
        init_pair(2, COLOR_WHITE, COLOR_CYAN);
    }
#endif
    
    // allocate the grid
    grid = AllocateGrid();
    if(grid == NULL)
    {
        printf("ERROR: Could not allocate memory for grid\n");
        return -1;
    }
    
    // initialize and print the grid
    InitGrid(grid);
#ifdef PRINT_GRID
    PrintGrid(grid);
#endif
    
    for(int time = 1; time <= TIME_STEPS; ++time)
    {
        // sleep to observe updates
#ifdef PRINT_GRID
        usleep(SLEEP_TIME);
#endif
        
        // allocate temporary grid for storing the simultaneous update
        int** tempGrid = AllocateGrid();
        if(tempGrid == NULL)
        {
            printf("ERROR: Could not allocate memory for grid\n");
            return -1;
        }
    
        for(int i = 0; i < GRID_SIZE; ++i)
        {
            for(int j = 0; j < GRID_SIZE; ++j)
            {
                int numNeighborAlive = 0;
            
                for(int k = -1; k <= 1; ++k)
                {
                    for(int l = -1; l <= 1; ++l)
                    {
                        int row = (i + k) % GRID_SIZE;
                        int col = (j + l) % GRID_SIZE;
                        
                        if(row == -1)
                            row = GRID_SIZE - 1;
                        if(col == -1)
                            col = GRID_SIZE - 1;
                    
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
        
        free(grid);
        grid = tempGrid;
        
#ifdef USE_CURSES
        move(0, 0);
#endif
#ifdef PRINT_GRID
    PrintGrid(grid);
#endif
    }
    
    CountLive(grid);
    
#ifdef USE_CURSES
    printw("\nPress any key to continue...");
    refresh();
    getchar();
#endif
    
    free(grid);
#ifdef USE_CURSES
    endwin();
#endif

    return 0;
}

int** AllocateGrid()
{
    int** grid;
    
    grid = malloc(sizeof(int) * GRID_SIZE);
    
    if(grid == NULL)
    {
        printf("ERROR: Could not allocate memory for grid\n");
        return NULL;
    }
    
    for(int i = 0; i < GRID_SIZE; ++i)
    {
        grid[i] = malloc(sizeof(int) * GRID_SIZE);
        
        if(grid == NULL)
        {
            printf("ERROR: Could not allocate memory for grid\n");
            return NULL;
        }
    }
    
    return grid;
}

void InitGrid(int** grid)
{
    for(int i = 0; i < GRID_SIZE; ++i)
    {
        for(int j = 0; j < GRID_SIZE; ++j)
        {
            grid[i][j] = rand() % 2;
        }
    }
}

void PrintGrid(int** grid)
{
    for(int i = 0; i < GRID_SIZE; ++i)
    {
        for(int j = 0; j < GRID_SIZE; ++j)
        {
#ifdef USE_CURSES
            if(has_colors())
            {
                addch(grid[i][j] == 1 ? (' ' | COLOR_PAIR(2)) : (' ' | COLOR_PAIR(1)));
                addch(grid[i][j] == 1 ? (' ' | COLOR_PAIR(2)) : (' ' | COLOR_PAIR(1)));
            }
            else
                addch(grid[i][j] == 1 ? '1' : '0');
            
            if(j == GRID_SIZE-1)
                addch('\n');
#else
            putchar(grid[i][j] == 1 ? '1' : '0');
            
            if(j == GRID_SIZE-1)
                putchar('\n');
#endif
        }
    }

#ifdef USE_CURSES
    refresh();
#else
    putchar('\n');
#endif
}

void CountLive(int** grid)
{
    int numLive = 0;

    for(int i = 0; i < GRID_SIZE; ++i)
    {
        for(int j = 0; j < GRID_SIZE; ++j)
        {
            numLive += grid[i][j];
        }
    }
    
#ifdef USE_CURSES
    printw("\nNumber of live cells = %d\n", numLive);
    refresh();
#else
    printf("Number of live cells = %d\n", numLive);
#endif
}

