#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>
#include <sys/time.h>
#include <malloc.h>

const int HOSTNAME_WITH_TIME_MESSAGE_TAG= 42;
const int HOSTNAME_LENGTH = 128;
const int MASTER_PROCESS = 0;

struct hostname_with_time {
  struct timeval tv;
  char hostname[128];
};

int main(int argc, char **argv)
{
	//Initialisiere die MPI-Umgebung mit den Parametern der Funktion
  MPI_Init(&argc, &argv);

  //Let's see how many Processes we got
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  //And what Process we are
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  //Single Process? Boring!
  if (world_size < 2) {
    fprintf(stderr,"Requires at least two processes.\n");
    exit(-1);
  }

  //The MPI Datatypes to declare
  MPI_Datatype MPI_HOSTNAME_WITH_TIME, MPI_HOSTNAME, MPI_TIMEVAL;

  //Hostname Type
  MPI_Type_contiguous(HOSTNAME_LENGTH, MPI_CHAR, &MPI_HOSTNAME);

  //Timeval Type
  int timeval_blocklengths[2] = {4,4};
  MPI_Datatype timeval_subtypes[2] = {MPI_LONG, MPI_LONG};
  MPI_Aint timeval_offsets[2];

  timeval_offsets[0] = offsetof(struct timeval, tv_sec);
  timeval_offsets[1] = offsetof(struct timeval, tv_usec);

  MPI_Type_struct(1, timeval_blocklengths, timeval_offsets, timeval_subtypes, &MPI_TIMEVAL);
  MPI_Type_commit(&MPI_TIMEVAL);

  //Hostname with Time Type
  int hostname_with_time_blocklengths[2] = {8, 128};
  MPI_Datatype hostname_with_time_subtypes[2] = {MPI_TIMEVAL, MPI_HOSTNAME};
  MPI_Aint hostname_with_time_offsets[2];

  hostname_with_time_offsets[0] = offsetof(struct hostname_with_time, tv);
  hostname_with_time_offsets[1] = offsetof(struct hostname_with_time, hostname);

  MPI_Type_struct(1, hostname_with_time_blocklengths, hostname_with_time_offsets, hostname_with_time_subtypes, &MPI_HOSTNAME_WITH_TIME);
  MPI_Type_commit(&MPI_HOSTNAME_WITH_TIME);

  if (world_rank != MASTER_PROCESS)
  {
    sleep(1);

    struct hostname_with_time *p = (struct hostname_with_time*) malloc(sizeof (struct hostname_with_time));

    gethostname(p->hostname, sizeof(p->hostname));

    gettimeofday(&(p->tv), NULL);

    printf("%d Sending now\n",world_rank);
    MPI_Send((void *) &p, 1, MPI_HOSTNAME_WITH_TIME, MASTER_PROCESS, HOSTNAME_WITH_TIME_MESSAGE_TAG, MPI_COMM_WORLD);
  }

 	//Prozess 0 gibt den String aus
  	//Nach Rang der Prozesse geordnet
  if (world_rank == MASTER_PROCESS)
  {
   for (int i = 1; i < world_size; ++i)
   {
    void *buf = malloc(sizeof(struct hostname_with_time));
    printf("%d Waiting for %d\n", world_rank,i);
    //MPI_Status *status;
    MPI_Recv(buf, 1, MPI_HOSTNAME_WITH_TIME, i, HOSTNAME_WITH_TIME_MESSAGE_TAG, MPI_COMM_WORLD, NULL);
    printf("%d Bla2\n", world_rank);
    struct hostname_with_time *message = (struct hostname_with_time*) buf;
    printf("%s\n", message->hostname);
  }
}

  //Prozesse beenden erst wenn alle Ausgaben fertig sind
printf("Rang %d waiting at Barrier\n", world_rank);
  MPI_Barrier(MPI_COMM_WORLD);
printf("Rang %d hopped over Barrier\n", world_rank);


  	//Rang X beendet jetzt!
printf("Rang %d beendet jetzt\n", world_rank);

  	//Beende die MPI-Umgebung
    MPI_Finalize();

    return 0;
}