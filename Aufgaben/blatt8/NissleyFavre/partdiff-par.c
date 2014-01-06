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
#include <mpi.h>

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

//changed: variables for the parallelization
//start and end of the matrixpart of the process in the whole matrix
int start;
int end;
int my_rank;
int number_of_processes;
int lines_in_mymatrix;

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

	//changed: calculate start and end
	int rest = (arguments->N + 1) % number_of_processes;
	//give every process the same number of lines
	lines_in_mymatrix = (arguments->N + 1)/number_of_processes;
	//distribute the rest over the first processes
	int rest_offset = 0;

	// compute offset of all smaller ranks
	for (int i = 0; i < my_rank; ++i)
	{
		if (i < rest)
		{
			rest_offset++;
		}
	}
	// compute start rank + all offsets off the smaller ranks 
	start = my_rank * lines_in_mymatrix + rest_offset;

	//offset of my rank
	if(my_rank < rest) 
	{
		rest_offset++;
	}

	end = (my_rank + 1) * lines_in_mymatrix + rest_offset - 1;

	//give every process his first and his last line, that he doesnt compute
	//but needs from the other processes to compute his lines
	if (my_rank > 0)
	{
		start--;
		lines_in_mymatrix++;
	}

	if (my_rank < (number_of_processes - 1))
	{
		end++;
		lines_in_mymatrix++;
	}

	if (my_rank < rest)
	{
		lines_in_mymatrix++;
	}

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

	arguments->M = allocateMemory(arguments->num_matrices * (lines_in_mymatrix) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory(lines_in_mymatrix * sizeof(double*));

		for (j = 0; j < lines_in_mymatrix; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (lines_in_mymatrix)) + (j * (N + 1));
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
		for (i = 0; i < lines_in_mymatrix; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}
	/* initialize borders, depending on function (function 2: nothing to do) */
	//changed: initialise the borders of your part of the matrix
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
			for (i = 0; i <= N; i++)
			{
				//initialise the first line
				if (my_rank == 0)
				{
					Matrix[g][0][i] = 1.0 - (h * i);
				}
				//initialise the last line
				if (my_rank == (number_of_processes - 1))
				{
					Matrix[g][lines_in_mymatrix - 1][i] = h * i;
					Matrix[g][lines_in_mymatrix - 1][0] = 0.0;
				}
			}
			for (i = 0; i < lines_in_mymatrix; ++i)
			{	
				//initialise all the lines between the first and the last
				Matrix[g][i][0] = 1.0 - (h * (i + start));
				Matrix[g][i][N] = h * (i + start);

			}

			Matrix[g][lines_in_mymatrix - 1][0] = 0.0;
			Matrix[g][0][N] = 0.0;
		}
	}
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate_jacobi (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options)
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

	//changed: build an array for the requests
	MPI_Request requests_before[2];
	MPI_Request requests_next[2];

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

		//changed: Send and Receive Messages

		//calculate the second line and send it to the process -1
		double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)(1 + start));
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[0][j] + Matrix_In[1][j-1] + Matrix_In[1][j+1] + Matrix_In[2][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[1][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[1][j] = star;
			}
		if (my_rank > 0)
		{
			//send the second line to the process - 1
			MPI_Isend(&Matrix_Out[1][0], N + 1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, &requests_before[0]);
		}

		//calculate the second last line and send it to the next process
		fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)(lines_in_mymatrix - 2 + start));
			}

			/* over all columns */
			for (j = 1; j < N; j++)
			{
				star = 0.25 * (Matrix_In[lines_in_mymatrix - 3][j] + Matrix_In[lines_in_mymatrix - 2][j-1] + Matrix_In[lines_in_mymatrix - 2][j+1] + Matrix_In[lines_in_mymatrix - 1][j]);

				if (options->inf_func == FUNC_FPISIN)
				{
					star += fpisin_i * sin(pih * (double)j);
				}

				if (options->termination == TERM_PREC || term_iteration == 1)
				{
					residuum = Matrix_In[lines_in_mymatrix - 2][j] - star;
					residuum = (residuum < 0) ? -residuum : residuum;
					maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
				}

				Matrix_Out[lines_in_mymatrix - 2][j] = star;
			}

		if (my_rank != number_of_processes - 1)
		{
			//send the second last line to the next process
			MPI_Isend(Matrix_Out[lines_in_mymatrix - 2], N + 1, MPI_DOUBLE, my_rank + 1, 1, MPI_COMM_WORLD, &requests_next[0]);
		}
		//be ready to receive the first and the last line
		if (my_rank > 0)
		{
			MPI_Irecv(Matrix_Out[0], N + 1, MPI_DOUBLE, my_rank - 1, MPI_ANY_TAG, MPI_COMM_WORLD, &requests_before[1]);
		}
		if (my_rank != number_of_processes - 1)
		{
			MPI_Irecv(Matrix_Out[lines_in_mymatrix - 1], N + 1, MPI_DOUBLE, my_rank + 1, MPI_ANY_TAG, MPI_COMM_WORLD, &requests_next[1]);
		}
		/* over all rows */
		//changed: do the calculation only for your part of the matrix
		for (i = 2; i < lines_in_mymatrix - 2; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)(i + start));
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
		MPI_Status s;

		if (my_rank > 0)
		{
			//wait until all messages to and from process -1 are received
			MPI_Waitall(2, requests_before, &s);
		}

		if (my_rank < number_of_processes - 1)
		{
			//wait until all messages to and from the next process are received
			MPI_Waitall(2, requests_next, &s);
		}
	
		/* check for stopping calculation, depending on termination method */
		if (options->termination == TERM_PREC)
		{
			//changed: calculate the maximum of all local maxresiduums
			MPI_Allreduce(&maxresiduum, &maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			if (maxresiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	results->m = m2;

	MPI_Allreduce(&maxresiduum, &maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	results->stat_precision = maxresiduum;
	
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void 
calculate_gauss (struct calculation_arguments const* arguments, struct calculation_results *results, struct options const* options)
{
	int i, j;                                   /* local variables for loops  */
	int m1, m2;                                 /* used as indices for old and new matrices       */
	double star;                                /* four times center value minus 4 neigh.b values */
	double residuum;                            /* residuum of current iteration                  */
	double maxresiduum;                         /* maximum residuum value of a slave in iteration */
	//changed: global maxresiduum
	double maxresiduum_global;
	//changed: variable for the abort signal
	int message_received;
	int abort_tag = 4;

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	/* initialize m1 and m2 depending on algorithm */
	//changed: only do what you need for gauß seidel
	m1 = 0;
	m2 = 0;

	//changed: build an array for the requests
	MPI_Request requests_before[2];
	MPI_Request requests_next[2];

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	while (term_iteration > 0)
	{
		double** Matrix_Out = arguments->Matrix[m1];
		double** Matrix_In  = arguments->Matrix[m2];
		MPI_Status s;
		maxresiduum = 0;
		maxresiduum_global = 0;

		if (options->termination == TERM_PREC)
		{
			//check if there is an abort message
			if (my_rank == 0)
			{
				MPI_Iprobe(number_of_processes - 1, abort_tag, MPI_COMM_WORLD, &message_received, &s);
			} 
			else 
			{
				MPI_Iprobe(my_rank -1, abort_tag, MPI_COMM_WORLD, &message_received, &s);
			}
			//changed: if you receive an abort-message: 
			if (message_received == 1)
			{
				printf("abort-message received from %d\n", my_rank);

				if (my_rank == 0)
				{
					MPI_Recv(&abort_tag, 1, MPI_INT, number_of_processes - 1, abort_tag, MPI_COMM_WORLD, &s);
				}
				else
				{
					MPI_Recv(&abort_tag, 1, MPI_INT, my_rank - 1, abort_tag, MPI_COMM_WORLD, &s);
				}

				term_iteration = 0;

				if (my_rank < number_of_processes -1)
				{
					MPI_Send(&abort_tag, 1, MPI_INT, my_rank + 1, abort_tag, MPI_COMM_WORLD);
					//MPI_Recv(Matrix_Out[lines_in_mymatrix - 1], N + 1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD, &s);
				}
			}
		}

		//changed: receive the first line and the maxresiduum
		if (my_rank > 0)
		{
			MPI_Recv(Matrix_Out[0], N + 1, MPI_DOUBLE, my_rank - 1, 1, MPI_COMM_WORLD, NULL);

			if (term_iteration == 0)
			{
				printf("received my first line from rank %d\n", my_rank);
			}

			MPI_Irecv(&maxresiduum_global, 1, MPI_DOUBLE, my_rank - 1, 2, MPI_COMM_WORLD, &requests_before[1]);
		}

		//changed: calculate the second line and send it to the process - 1
		double fpisin_i = 0.0;

		if (options->inf_func == FUNC_FPISIN)
		{
			fpisin_i = fpisin * sin(pih * (double)(1 + start));
		}

		/* over all columns */
		for (j = 1; j < N; j++)
		{
			star = 0.25 * (Matrix_In[0][j] + Matrix_In[1][j-1] + Matrix_In[1][j+1] + Matrix_In[2][j]);

			if (options->inf_func == FUNC_FPISIN)
			{
				star += fpisin_i * sin(pih * (double)j);
			}

			if (options->termination == TERM_PREC || term_iteration == 1)
			{
				residuum = Matrix_In[1][j] - star;
				residuum = (residuum < 0) ? -residuum : residuum;
				maxresiduum = (residuum < maxresiduum) ? maxresiduum : residuum;
			}

			Matrix_Out[1][j] = star;
		}

		if (my_rank > 0 && term_iteration != 0)
		{
			//send the second line to the process - 1
			MPI_Isend(&Matrix_Out[1][0], N + 1, MPI_DOUBLE, my_rank - 1, 0, MPI_COMM_WORLD, &requests_before[0]);
		}

		/* over all rows */
		//calculate from the third to the second last line in your part of the matrix
		for (i = 2; i < lines_in_mymatrix - 1; i++)
		{
			double fpisin_i = 0.0;

			if (options->inf_func == FUNC_FPISIN)
			{
				fpisin_i = fpisin * sin(pih * (double)(i + start));
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

		//changed: send the second last line to the next process
		if (my_rank < number_of_processes - 1)
		{
			MPI_Isend(Matrix_Out[lines_in_mymatrix - 2], N + 1, MPI_DOUBLE, my_rank + 1, 1, MPI_COMM_WORLD, &requests_next[0]);
			if (term_iteration == 0)
			{
				printf("send the next process his first line\n");
			}
		}

		results->stat_iteration++;
		results->stat_precision = maxresiduum;

		/* exchange m1 and m2 */
		i = m1;
		m1 = m2;
		m2 = i;

		if (my_rank > 0 && term_iteration != 0)
		{
			//wait until process - 1 has received his last line and I have received the maxresiduum
			MPI_Waitall(2, requests_before, &s);
		}

		if (my_rank < number_of_processes - 1)
		{
			//wait until the next process has received his first line
			MPI_Waitall(1, requests_next, &s);
		}

		//receive my last line from the next process
		if (my_rank != number_of_processes - 1 && term_iteration != 0)
		{
			MPI_Recv(Matrix_Out[lines_in_mymatrix - 1], N + 1, MPI_DOUBLE, my_rank + 1, 0, MPI_COMM_WORLD, &s);
		}

		//reduce the maxresiduum
		maxresiduum_global = (maxresiduum < maxresiduum_global) ? maxresiduum_global : maxresiduum;
		if (my_rank < number_of_processes - 1)
		{
			MPI_Send(&maxresiduum_global, 1, MPI_DOUBLE, my_rank + 1, 2, MPI_COMM_WORLD);
		}

		//check for termination depending on the termination-method 
		if (options->termination == TERM_PREC)
		{
			if (my_rank == number_of_processes - 1 && term_iteration != 0)
			{
				if (maxresiduum_global < options->term_precision)
				{
					MPI_Send(&abort_tag, 1, MPI_INT, 0, abort_tag, MPI_COMM_WORLD);
				}
			}


		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}

	MPI_Allreduce(&maxresiduum, &maxresiduum, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
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

void DisplayMatrix ( char *s, double *v, int interlines , int rank , int size, int from, int to )
{
  int x,y;
  int lines = 8 * interlines + 9;
  MPI_Status status;

  /* first line belongs to rank 0 */
  if (rank == 0)
    from--;

  /* last line belongs to rank size-1 */
  if (rank + 1 == size)
    to++;

  if (rank == 0)
    printf ( "%s\n", s );

  for ( y = 0; y < 9; y++ )
  {
    int line = y*(interlines+1);

    if (rank == 0)
    {
      /* check whether this line belongs to rank 0 */
      if (line < from || line > to)
      {
        /* use the tag to receive the lines in the correct order
         * the line is stored in v, because we do not need it anymore */
        MPI_Recv(v, lines, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
      }
    }
    else
    {
      if (line >= from && line <= to)
      {
        /* if the line belongs to this process, send it to rank 0
         * (line - from + 1) is used to calculate the correct local address */
        MPI_Send(&v[(line - from + 1)*lines], lines, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
      }
    }

    for ( x = 0; x < 9; x++ )
    {
      if (rank == 0)
      {
        if (line >= from && line <= to)
        {
          /* this line belongs to rank 0 */
          printf ( "%7.4f", v[line*lines + x*(interlines+1)]);
        }
        else
        {
          /* this line belongs to another rank and was received above */
          printf ( "%7.4f", v[x*(interlines+1)]);
        }
      }
    }

    if (rank == 0)
      printf ( "\n" );
  }
  fflush ( stdout );
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	//changed: initialise the MPI-Environment
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	/* get parameters */
	AskParams(&options, argc, argv);              /* ************************* */

	initVariables(&arguments, &results, &options);           /* ******************************************* */
	
	allocateMatrices(&arguments);        /*  get and initialize variables and matrices  */

	initMatrices(&arguments, &options);            /* ******************************************* */

	gettimeofday(&start_time, NULL);                   /*  start timer         */
	if (options.method == METH_JACOBI)
	{
		calculate_jacobi(&arguments, &results, &options);
	}   
	else if (options.method == METH_GAUSS_SEIDEL)
	{
		calculate_gauss(&arguments, &results, &options);
	}                              /*  solve the equation  */
	gettimeofday(&comp_time, NULL);                   /*  stop timer          */

	displayStatistics(&arguments, &results, &options);
	DisplayMatrix("Result Matrix: \n", &arguments.Matrix[results.m][0][0], options.interlines , my_rank, number_of_processes, start + 1, end - 1);

	freeMatrices(&arguments);                                       /*  free memory     */

	MPI_Finalize();

	return 0;
}
