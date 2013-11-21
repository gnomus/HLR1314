#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <time.h>


int main(int argc, char const **argv)
{
	//Initialisiere die MPI-Umgebung mit den Parametern der Funktion
	MPI_Init(&argc, &argv);

	//Anzahl der Prozesse
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	//Rang des Prozesses
  	int world_rank;
  	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  	//Prozess 1-n erzeugen HOSTNAME: TIMESTAMP als string und 
  	//übergeben diesen an Prozess 0
  	if (world_rank != 0)
  	{
  		char *hostname;
  		gethostname(hostname, sizeof(hostname));
 		char *timestamp; 
 		gettimeofday(timestamp, NULL);
 		char *ausgabe = hostname + ":" + timestamp;

 		MPI_Send((void *) ausgabe, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD)
	}

 	//Prozess 0 gibt den String aus
  	//Nach Rang der Prozesse geordnet
 	if (world_rank == 0)
 	{
 		for (int i = 1; i < world_size; ++i)
 		{
			MPI_Recv(ausgabe, 1, MPI_CHAR, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
 			printf("%s:%s\n",(char *) ausgabe);
 		}
	}

  	//Prozesse beenden erst wenn alle Ausgaben fertig sind
  	MPI_Barrier(MPI_COMM_WORLD);

  	//Rang X beendet jetzt!
  	printf("Rang %d beendet jetzt\n", world_rank);

  	//Beende die MPI-Umgebung
  	MPI_Finalize();

	return 0;
}