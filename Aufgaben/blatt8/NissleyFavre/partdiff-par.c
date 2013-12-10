/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                TU Muenchen - Institut fuer Informatik                  **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**            JK und andere  besseres Timing, FLOP Berechnung             **/
/**                                                                        **/
/** File:      partdiff-seq.c                                              **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauss-Seidel and   **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>

#include "partdiff-par.h"

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */

//changed: Variables for Parallelization
int myrank;
int number_of_ranks; 
int size_of_mymatrix; 
int start;
int end; 

/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	//Anzahl der Zeilen, auch gebraucht von jedem Prozess
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	//changed: berechne size_of_mymatrix
	calculate_size_of_my_matrix(arguments->N);

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}

/* ************************************************************************ */
/* calculate the size of the matrix depending on the rank of the process    */
/* ************************************************************************ */
static
void
calculate_size_of_my_matrix(uint64_t size_of_whole_matrix)
{
	//give every process the same number of lines
	lines_per_process = N/number_of_ranks;
	rest = N % number_of_ranks;
	size_of_mymatrix = lines_per_process;

	int restezusatz = 0;
	//distribute the rest
	if (myrank < rest)
	{
		size_of_mymatrix++;
		start = myrank * size_of_mymatrix + restezusatz;
		end = (myrank + 1) * size_of_mymatrix + restezusatz - 1;
		restezusatz++;
	}

	if (myrank > rest)
	{
		start = myrank * size_of_mymatrix + restezusatz;
		end = (myrank + 1) * size_of_mymatrix + restezusatz - 1;
	}

	//give every process two extra lines for data sharing
	size_of_matrix_0 = size_of_mymatrix + 1;
	size_of_mymatrix = size_of_mymatrix + 2;

	if (myrank == 0 || myrank == N - 1)
	{
		size_of_mymatrix = size_of_matrix_0;
	}
}


/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
static
void
freeMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes)\n", size);
		/* exit program */
		exit(1);
	}

	return p;
}

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;

	//changed: only allocate the space one process needs
	arguments->M = allocateMemory(arguments->num_matrices * size_of_mymatrix * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory(size_of_mymatrix * sizeof(double*));

		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
	
	//at first, initialise all fields with 0 
	uint64_t g, i, j;                                /*  local variables for loops   */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= size_of_mymatrix; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}


	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
			//changed: initialize borders depending on rank of the process
			if (myrank == 0)
			{
				Matrix[g][0][N] = 0.0;

				for (int i = 1; i < N; ++i)
				{
					//initalize line 0
					Matrix[g][0][i] = 1.0 - (h * i);

					//todo: zweites i geht nur bis sizeofmymatrix und beginnt bei 1
					//initialize first and last element of the lines between 0 and N 
					//das i muss auf die Gesamtmatrix bezogen sein
					Matrix[g][i][0] = 1.0 - (h * i);
					Matrix[g][i][N] = h * i;
				}
			}

			if (myrank != 0 && myrank != number_of_ranks)
			{
				
				for (int i = 0; i < N; ++i)
				{
					//initialize first and last element of the lines between 0 and N 
					//todo
					Matrix[g][i][0] = 1.0 - (h * (i + start);
					Matrix[g][i][N] = h * (i + start);
				}
			}

			if (myrank == (number_of_ranks - 1))
			{
				Matrix[g][size_of_mymatrix - 1][0] = 0.0;

				for (int i = 0; i < (N-1); ++i)
				{
					//initialize last line
					Matrix[g][size_of_mymatrix - 1][i] = h * (i + start);

					//initialize first and last element of the lines between 0 and N 
					//todo
					Matrix[g][i][0] = 1.0 - (h * (i + start));
					Matrix[g][i][N] = h * (i + start);
				}
			}
		}
	}
	
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options)
{
	int i, j;                                   /* local variables for loops  */
	int m1, m2;                                 /* used as indices for old and new matrices       */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration                  */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];

		maxresiduum = 0;

		/* over all rows */
		//changed: Berechnung nur in der Teilmatrix, erste und letzte Zeile ausgelassen
		for (i = 1; i < size_of_mymatrix; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[i][j] = star;
			}
		}

		results->stat_iteration++;
		results->stat_precision = maxresiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation, depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxresiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}

		//changed: the process needs to wait until all other processes have finished
		 MPI_Barrier(MPI_COMM_WORLD);

		 if (myrank % 2 == 0)
		 {
		 	if (myrank != 0)
		 	{
		 		//Empfängt erste Zeile von Prozess i - 1
		 		MPI_Recv (void* buf, int count, MPI_Datatype datatype, myrank - 1, int tag, MPI_Comm comm, MPI_Status* status);
		 		//Sendet erste zeile nach Prozess i - 1
		 		MPI_Send(void* buf, int count, MPI_Datatype datatype, myrank - 1, int tag, MPI_COMM_WORLD);
		 	}

		 	if (myrank != number_of_ranks)
		 	{
		 		//Sendet letzte zeile nach Prozess i + 1
		 		MPI_Send(void* buf, int count, MPI_Datatype datatype, myrank + 1, int tag, MPI_COMM_WORLD);
		 		//Empfängt erste Zeile von Prozess i + 1
		 		MPI_Recv (void* buf, int count, MPI_Datatype datatype, myrank + 1, int tag, MPI_Comm comm, MPI_Status* status);
		 	}
		 }

		 //todo: Einbau der empfangenen Zeilen
	}

	results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauss-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y) = 0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC)
	{
		printf("Hinreichende Genaugkeit");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion DisplayMatrix:                               **/
/**                                                                        **/
/** Die Funktion DisplayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			printf ("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
		}

		printf ("\n");
	}

	fflush (stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	//changed: Start der Parallelisierung
	MPI_Init(&argc, &argv);

	//changed: initialize variables for parallelisation
	MPI_Comm_size(MPI_COMM_WORLD, &number_of_ranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	/* get parameters */
	AskParams(&options, argc, argv);              /* ************************* */

	initVariables(&arguments, &results, &options);           /* ******************************************* */

	allocateMatrices(&arguments);        /*  get and initialize variables and matrices  */
	initMatrices(&arguments, &options);            /* ******************************************* */

	gettimeofday(&start_time, NULL);                   /*  start timer         */
	calculate(&arguments, &results, &options);                                      /*  solve the equation  */
	gettimeofday(&chomp_time, NULL);                   /*  stop timer          */

	displayStatistics(&arguments, &results, &options);
	//changed: use the display function for the mpi programm
	// todo: to und from für jeden Prozess
	DisplayMatrix( "Matrix:\n", calculation_arguments->M, options->interlines, myrank, number_of_ranks, start, end);

	freeMatrices(&arguments);                                       /*  free memory     */

	//changed: Ende der Parallelisierung
	MPI_Finalize();

	return 0;
}
