#include <stdio.h>
#include <mpi.h>


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

  	//Name des Prozessors
  	char processor_name[MPI_MAX_PROCESSOR_NAME];
  	int name_len;
  	MPI_Get_processor_name(processor_name, &name_len);


  	//Prozess 1-n erzeugen HOSTNAME: TIMESTAMP als string und 
  	//übergeben diesen an Prozess 0
  	if (world_rank != 0)
  	{
  		char *hostname = gethostname();
 		char *timestamp = gettimeofday();
 		char *ausgabe = hostname + ":" + timestamp;

 		MPI_Send((void *) ausgabe, int count, MPI_Datatype datatype, 0, int tag, MPI_Comm communicator)
	}
  	
 	//Prozess 0 gibt den String aus
  	//Nach Rang der Prozesse geordnet
 	if (world_rank == 0)
 	{
 		for (int i = 0; i < world_size; ++i)
 		{
			MPI_Recv(ausgabe, int count, MPI_Datatype datatype, i, int tag, MPI_Comm communicator, MPI_Status* status);
 			printf("%s:%s\n",(char *) ausgabe);
 		}
	}

  	//Prozesse beenden erst wenn alle Ausgaben fertig sind
  	
  	MPI Barrier();

  	//Rang X beendet jetzt!
  	printf("Rang %d beendet jetzt\n", world_rank);

  	//Beende die MPI-Umgebung
  	MPI_Finalize();

	return 0;
}