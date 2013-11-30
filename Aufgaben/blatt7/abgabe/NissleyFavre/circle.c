#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define MASTER_PROCESS 0

int* init (int N, int bufsize, int count, int rest)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  printf("rank %i | N %i | bufsize %i | count %i | rest %i\n", rank, N, bufsize, count, rest);

  int *buf = (int *) malloc(sizeof(int) * bufsize);

  srand(time(NULL));

  for (int i = 0; i < count; i++)
  {
    buf[i] = rand() % 25; //do not modify %25
    printf("%d\n", buf[i]);
  }

  if(rank >= rest) {
    buf[bufsize] = -1;
  }

  return buf;
}

int* circle (int* buf)
{
  //todo
  return buf;
}

int
main (int argc, char** argv)
{
  unsigned *uN;
  int rank, world_size;
  int* buf;

    //CHANGED: Initialize MPI
  MPI_Init(&argc, &argv);
  //CHANGED: Rank
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if (argc < 2)
  {
    if (rank == MASTER_PROCESS) {
      printf("Arguments error\n");
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }


  if (sscanf(argv[1], "%u", uN) != 1) {
    if (rank == MASTER_PROCESS) {
      printf("Please specify Array length\n");
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  int N = *uN;

  //Abbruch bei zu vielen Prozessen
  if(world_size > N)
  {
    if (rank == MASTER_PROCESS)
    {
      printf("Error: More Processes than Array Elements\n");
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  int count = N / world_size;
  int rest = N % world_size;

  int bufsize = (rest != 0) ? count + 1 : count;
  if (rank < rest)
  {
    count++;
  }

  buf = init(N, bufsize, count, rest);


  printf("\nBEFORE\n");

  for (int i = 0; i < bufsize; i++)
  {
    //printf ("rank %i: %i\n", rank, buf[i]);
  }

  circle(buf);

  //printf("\nAFTER\n");

  for (int j = 0; j < bufsize; j++)
  {
    //printf ("rank %d: %d\n", rank, buf[j]);
  }


  MPI_Finalize();
  return EXIT_SUCCESS;
}
