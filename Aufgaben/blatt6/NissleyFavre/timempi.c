#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
#include <malloc.h>

const int HOSTNAME_WITH_TIME_MESSAGE_TAG = 42;
const int MASTER_PROCESS = 0;
const int HOSTNAME_WITH_TIME_MESSAGE_SIZE = 256;


int main(int argc, char **argv)
{
  //Initializing MPI
  MPI_Init(&argc, &argv);

  //Let's see how many Processes we got
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  //And what Process we are
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  //Single Process? Boring!
  if (world_size < 2) {
    if (world_rank == MASTER_PROCESS)
      fprintf(stderr,"Requires at least two processes.\n");
    MPI_Finalize();
    exit(-1);
  }

  //We are not the Master Process. Let's do this.
  if (world_rank != MASTER_PROCESS)
  {
    //Gathering Hostname
    char hostname[128];
    gethostname(hostname, sizeof(hostname));

    //What Time is it?
    struct timeval tv;
    time_t nowtime;
    struct tm *nowtm;
    char tmbuf[64], tbuf[64];

    gettimeofday(&tv, NULL);
    nowtime = tv.tv_sec;
    nowtm = localtime(&nowtime);

    //Stringbuilding in C is fun!
    strftime(tmbuf, sizeof tmbuf, "%Y-%m-%d %H:%M:%S", nowtm);
    snprintf(tbuf, sizeof tbuf, "%s.%06li", tmbuf, tv.tv_usec);

    //More Stringbuilding..
    char buf[HOSTNAME_WITH_TIME_MESSAGE_SIZE];
    snprintf(buf, sizeof buf, "%s: %s", hostname, tbuf);

    //Let's send this to the Master Process!
    MPI_Send((void *) buf, HOSTNAME_WITH_TIME_MESSAGE_SIZE, MPI_CHAR, MASTER_PROCESS, HOSTNAME_WITH_TIME_MESSAGE_TAG, MPI_COMM_WORLD);

  }

  //We are the master Process. Let's do something else
  if (world_rank == MASTER_PROCESS)
  {
    //We need some ressources for saving recieved Messages in Memmory
    void *buf = malloc(sizeof(char)*HOSTNAME_WITH_TIME_MESSAGE_SIZE);

    //We want to reviede a Message from every other Process (1 - N)
    for (int i = 1; i < world_size; ++i)
    {
      //Recieving...
      MPI_Recv(buf, HOSTNAME_WITH_TIME_MESSAGE_SIZE, MPI_CHAR, i, HOSTNAME_WITH_TIME_MESSAGE_TAG, MPI_COMM_WORLD, NULL);
      //Printing...
      printf("%s\n",(char*)buf);
    }
    //Don't need those Ressources anymore
    free(buf);
  }

  //Let's be fair and wait for everyone :)
  MPI_Barrier(MPI_COMM_WORLD);
  //Notify the User about the soon to come exit
  printf("Rang %d beendet jetzt\n", world_rank);
  //MPI_Barrier(MPI_COMM_WORLD);

  //Close the MPI Environment
  MPI_Finalize();

  //Exit! Finally :D
  return 0;

}