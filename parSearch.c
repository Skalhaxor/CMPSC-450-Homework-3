// Parallel Search
// Food
//

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {
   if(argc <= 2) {
      fprintf(stderr, "Usage: %s filename target maxlines\n", argv[0]);
      return 1;
   }
  
   int err = MPI_Init(&argc, &argv);
   if(err != MPI_SUCCESS) {
      fprintf(stderr, "MPI_Init error\n");
      return 1;
   }   
   
   int rank, numProcs;
   int target;
   int dataLength;
   int* data;
   MPI_Status status;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
   
   // Read file section and send target value
   if(rank == 0) {
      // Some of this code should be error-checked, etc.
      // For example, malloc, fopen, and the MPI_ calls should be
      // But since this is just a rough practice, don't bother
      target = atoi(argv[2]);
      
      int maxLength = atoi(argv[3]);
      dataLength = 0;
      data = malloc(maxLength * sizeof(int));
      FILE* input = fopen(argv[1], "r");
      char line[12];
      while(dataLength  < maxLength && fgets(line, 11, input)) {
         data[dataLength] = atoi(line);
         dataLength++;
      }
   }

   MPI_Bcast(&target, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Bcast(&dataLength, 1, MPI_INT, 0, MPI_COMM_WORLD);
   printf("Processor %d of %d: %d!\n", rank, numProcs, target);

   // Send data segments
   int sendCounts[numProcs];
   int sendOffsets[numProcs];
   int segLength = dataLength / numProcs; 
   int extras = dataLength - segLength * numProcs;
   int recCount;
   for(int i = 0; i < numProcs; ++i) {
      sendCounts[i] = segLength;
      sendOffsets[i] = i * segLength;
   }
   sendCounts[numProcs - 1] += extras;
   recCount = sendCounts[rank];
   int dataSeg[recCount];
   MPI_Scatterv(data, sendCounts, sendOffsets, MPI_INT, dataSeg, recCount, MPI_INT, 0, MPI_COMM_WORLD);
   
   // Have each process search for its matches and print results
   for(int i=0; i < recCount; ++i) {
      if(dataSeg[i] == target) {
         printf("Process %d found a match at index %d\n", rank, sendOffsets[rank] + i);
      }
   }

   MPI_Finalize();

//   Non-parallel search
//   for(int i =0; i < length; ++i) {
//      if(data[i] == target) {
//         printf("Found on line %d\n", i+1);
//      }
//   }

}

