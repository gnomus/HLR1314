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
#include <pthread.h>

#include "partdiff-posix.h"


/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time;       /* time when program started                      */
struct timeval comp_time;        /* time when calculation completed                */
//changed: erstell mutex für das maxresiduum
pthread_mutex_t maxresiduum_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_barrier_t while_barrier;


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
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

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

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
	uint64_t g, i, j;                                /*  local variables for loops   */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
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
			for (i = 0; i <= N; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i);
				Matrix[g][i][N] = h * i;
				Matrix[g][0][i] = 1.0 - (h * i);
				Matrix[g][N][i] = h * i;
			}

			Matrix[g][N][0] = 0.0;
			Matrix[g][0][N] = 0.0;
		}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
//changed: calculate nimmt nurnoch eine struct entgegen, die alle parameter enthält weil den threads nur ein argument
//gegeben werden kann
static
void
calculate (void* opts)
{
	//Variablen die für jeden thread privat sind
	int i, j;                                   /* local variables for loops  */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration                  */
	int m1, m2;                                 /* used as indices for old and new matrices       */

    struct calculate_options* args  = (struct calculate_options*) opts;

	int const N = args->arguments->N;
	double const h = args->arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = args->options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	if (args->options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (args->options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = args->arguments->Matrix[m1];
		double** Matrix_In  = args->arguments->Matrix[m2];

		//changed: jeder thread bekommt einen start und einen endpunkt, so das hier jeder prozess
		//nur seinen teil abarbeitet, insgesamt wird 0 bis N bearbeitet
		/* over all rows */
		for (i = args->start; i <= args->ende; i++)
		{
			double fpisin_i = 0.0;

			if (args->options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)i);
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				/*printf("---------\n");
				printf("start: %i\n",args->start);
				printf("end: %i\n",args->ende);
				printf("N: %i\n",N);
				printf("i: %i\n",i);
				printf("j: %i\n",j);
				printf("---------\n");*/

				star = 0.25 * (Matrix_In[i-1][j] + Matrix_In[i][j-1] + Matrix_In[i][j+1] + Matrix_In[i+1][j]);

				if (args->options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (args->options->termination == TERM_PREC || term_iteration == 1)
				{
					pthread_mutex_lock(&maxresiduum_mutex);
					residuum = Matrix_In[i][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					args->maxresiduum = (residuum < args->maxresiduum) ? args->maxresiduum : residuum;
					pthread_mutex_unlock(&maxresiduum_mutex);
				}

				Matrix_Out[i][j] = star;
			}
		}

		//changed: globales maxresiduum aus lokalen maxresiduums berechnen
		//mutex damit die threads nicht gleichzeitig prüfen und schreiben
		//if (maxresiduum > args->maxresiduum)
		//{
		//	calculate_options->maxresiduum = maxresiduum;
		//}

		//Changed: Only let one Thread count the Iterations
		if (args->threadid == 0) args->results->stat_iteration++;

		args->results->stat_precision = args->maxresiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		/* check for stopping calculation, depending on termination method */
		if (args->options->termination == TERM_PREC)
		{
			if ( args->maxresiduum < args->options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (args->options->termination == TERM_ITER)
		{
			term_iteration--;
		}

		//changed: barrier eingefügt
		pthread_barrier_wait(&while_barrier);
	}

	args->results->m = m2;

	//changed: sorgt dafür dass die threads beendet werden
	pthread_exit(NULL);

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
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	/* get parameters */
	AskParams(&options, argc, argv);              /* ************************* */

	initVariables(&arguments, &results, &options);           /* ******************************************* */

	int N = arguments.N;

	pthread_barrier_init(&while_barrier, NULL, options.number);

	//changed: array für die threads erstellen
	pthread_t *threads = malloc(sizeof(pthread_t)*options.number);


	allocateMatrices(&arguments);        /*  get and initialize variables and matrices  */
	initMatrices(&arguments, &options);            /* ******************************************* */

	gettimeofday(&start_time, NULL);

	double maxresiduum = 0;
	double *maxresiduum_pointer = &maxresiduum;
	int sizeof_block = (int) (N-1)/options.number;
	int rest = 0;
	//changed: erstellen der threads, jeder thread bekommt als id die adresse von dem element des arrays in dem er steht
	//die gesamte funktion calculate wird aufgeteilt                  /*  start timer         */
	for(uint64_t i = 0; i < options.number; i++)
	{
		//changed: fasse alle argumente für die threads in einer struct zusammen weil nur ein
		//argument an die threads übergeben werden kann
		//es wird eine start und eine end-variable für jeden thread definiert
		struct calculate_options *calculate_options_thread = malloc(sizeof(struct calculate_options));

		calculate_options_thread->options = &options;
		calculate_options_thread->arguments = &arguments;
		calculate_options_thread->results = &results;
		calculate_options_thread->maxresiduum = *maxresiduum_pointer;   /* maximum residuum value of a slave in iteration */
		calculate_options_thread->threadid = i;

		calculate_options_thread->start = i*sizeof_block + rest + 1;

		if(i < ((N-1)%options.number))
		{
			calculate_options_thread->ende = (i+1)*sizeof_block + 1 + rest;
			rest = rest + 1;
		}
		else
		{
			calculate_options_thread->ende = (i+1)*sizeof_block + rest;
		}
		printf("i: %i\n", i);
		printf("Start: %i\n", calculate_options_thread->start);
		printf("End: %i\n", calculate_options_thread->ende);

		pthread_create(&threads[i], NULL, calculate, (void *) calculate_options_thread);
	}
	                                      /*  free memory     */
	//changed: joine alle threads wieder zu einem masterthread
	for(uint64_t i=0; i < options.number; i++)
    {
	  pthread_join(threads[i], NULL);
	}
	                                    /*  solve the equation  */
	gettimeofday(&comp_time, NULL);                   /*  stop timer          */

	displayStatistics(&arguments, &results, &options);
	DisplayMatrix(&arguments, &results, &options);

	freeMatrices(&arguments);
	pthread_barrier_destroy(&while_barrier);


	return 0;
}
