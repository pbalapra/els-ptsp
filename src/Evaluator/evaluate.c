/*
 
 ######################################################
 ########## Evaluator for the PTSP solutions ##########
 ######################################################
 
      Version: 1.0
      File:    evaluate.c
      Author:  Prasanna Balaprakash
      Purpose: Evaluating the PTSP solutions from PTSP algorithms.
      Note: The distance calculations are taken from ACOTSP code
      Copyright (C) 2008  Prasanna Balaprakash
*/


/***************************************************************************
 
    Program's name: Evaluator
 
    Evaluator for the results of the PTSP solutions 
 
	Copyright (C) 2008 Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Email: {pbalapra,mbiro,stuetzle}@ulb.ac.be
    Mail address:	Prasanna Balaprakash, 
					IRIDIA, Universite Libre de Bruxelles 
					50, Av. F. Roosevelt, CP 194/6 
					B-1050 Brussels, Belgium 
					http://iridia.ulb.ac.be/~prasanna

*****************************************************************************************/

#include <stdlib.h>
#include <argp.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <stddef.h>
#include <error.h>
#include "evaluate.h"
#include <gsl/gsl_math.h>


struct problem instance;


static double
dtrunc (double x)
{
  int k;

  k = (int) x;
  x = (double) k;
  return x;
}

long int
(*distance)(long int, long int);  /* function pointer */

/*
      FUNCTION: the following four functions implement different ways of 
                computing distances for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
*/

long int
round_distance (long int i, long int j)
/*
      FUNCTION: compute Euclidean distances between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
*/
{
  double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
  double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
  double r  = sqrt(xd*xd + yd*yd) + 0.5;

  return (long int) r;
}

long int ceil_distance (long int i, long int j)
/*
      FUNCTION: compute ceiling distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
*/
{
  double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
  double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
  double r  = sqrt(xd*xd + yd*yd) + 0.000000001;

  return (long int)r;
}

long int
geo_distance (long int i, long int j)
/*
      FUNCTION: compute geometric distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: adapted from concorde code
                for the definition of how to compute this distance see TSPLIB
*/
{
  double deg, min;
  double lati, latj, longi, longj;
  double q1, q2, q3;
  long int dd;
  double x1 = instance.nodeptr[i].x, x2 = instance.nodeptr[j].x;
  double y1 = instance.nodeptr[i].y, y2 = instance.nodeptr[j].y;

  deg = dtrunc (x1);
  min = x1 - deg;
  lati = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
  deg = dtrunc (x2);
  min = x2 - deg;
  latj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

  deg = dtrunc (y1);
  min = y1 - deg;
  longi = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
  deg = dtrunc (y2);
  min = y2 - deg;
  longj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

  q1 = cos (longi - longj);
  q2 = cos (lati - latj);
  q3 = cos (lati + latj);
  dd = (int) (6378.388 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
  return dd;

}

long int
att_distance (long int i, long int j)
/*
      FUNCTION: compute ATT distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
*/
{
  double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
  double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
  double rij = sqrt ((xd * xd + yd * yd) / 10.0);
  double tij = dtrunc (rij);
  long int dij;

  if (tij < rij)
    dij = (int) tij + 1;
  else
    dij = (int) tij;
  return dij;
}


struct point
      *read_ptsp(const char *ptsp_file_name)
      /*
            FUNCTION: parse and read instance file
            INPUT:    instance name
            OUTPUT:   list of coordinates for all nodes
            COMMENTS: Instance files have to be in PTSPLIB format, otherwise procedure fails
      */
  {
    FILE         *ptsp_file;
    char         buf[LINE_BUF_LEN];
    long int     i, j;
    struct point *nodeptr;


    ptsp_file = fopen(ptsp_file_name, "r");
    if ( ptsp_file == NULL )
      {
        fprintf(stderr,"No instance file specified, abort\n");
        exit(1);
      }
    assert(ptsp_file != NULL);
    /*printf("\nreading tsp-file %s ... \n\n", tsp_file_name);*/

    fscanf(ptsp_file,"%s", buf);
    while ( strcmp("NODE_COORD_SECTION", buf) != 0 )
      {
        if ( strcmp("NAME", buf) == 0 )
          {
            fscanf(ptsp_file, "%s", buf);
            TRACE ( printf("%s ", buf); )
            fscanf(ptsp_file, "%s", buf);
            strcpy(instance.name, buf);
            TRACE ( printf("%s \n", instance.name); )
            buf[0]=0;
          }
        else if ( strcmp("NAME:", buf) == 0 )
          {
            fscanf(ptsp_file, "%s", buf);
            strcpy(instance.name, buf);
            TRACE ( printf("%s \n", instance.name); )
            buf[0]=0;
          }
        else if ( strcmp("COMMENT", buf) == 0 )
          {
            fgets(buf, LINE_BUF_LEN, ptsp_file);
            TRACE ( printf("%s", buf); )
            buf[0]=0;
          }
        else if ( strcmp("COMMENT:", buf) == 0 )
          {
            fgets(buf, LINE_BUF_LEN, ptsp_file);
            TRACE ( printf("%s", buf); )
            buf[0]=0;
          }
        else if ( strcmp("TYPE", buf) == 0 )
          {
            fscanf(ptsp_file, "%s", buf);
            TRACE ( printf("%s ", buf); )
            fscanf(ptsp_file, "%s", buf);
            TRACE ( printf("%s\n", buf); )
            if( strcmp("PTSP", buf) != 0 )
              {
                fprintf(stderr,"\n Not a TSP or PTSP instance in TSPLIB format !!\n");
                exit(1);
              }
            buf[0]=0;
          }
        else if ( strcmp("TYPE:", buf) == 0 )
          {
            fscanf(ptsp_file, "%s", buf);
            TRACE ( printf("%s\n", buf); )
            if( strcmp("PTSP", buf) != 0 )
              {
                fprintf(stderr,"\n Not a TSP or PTSP instance in TSPLIB format !!\n");
                exit(1);
              }
            buf[0]=0;
          }
        else if( strcmp("DIMENSION", buf) == 0 )
          {
            fscanf(ptsp_file, "%s", buf);
            TRACE ( printf("%s ", buf); );
            fscanf(ptsp_file, "%ld", &n);
            instance.n = n;

            TRACE ( printf("%ld\n", n); );
            assert ( n > 2 && n < 20000);
            buf[0]=0;
          }
        else if ( strcmp("DIMENSION:", buf) == 0 )
          {
            fscanf(ptsp_file, "%ld", &n);
            instance.n = n;
            TRACE ( printf("%ld\n", n); );
            assert ( n > 2 && n < 6000);
            buf[0]=0;
          }
        else if( strcmp("DISPLAY_DATA_TYPE", buf) == 0 )
          {
            fgets(buf, LINE_BUF_LEN, ptsp_file);
            TRACE ( printf("%s", buf); );
            buf[0]=0;
          }
        else if ( strcmp("DISPLAY_DATA_TYPE:", buf) == 0 )
          {
            fgets(buf, LINE_BUF_LEN, ptsp_file);
            TRACE ( printf("%s", buf); );
            buf[0]=0;
          }
        else if( strcmp("EDGE_WEIGHT_TYPE", buf) == 0 )
          {
            buf[0]=0;
            fscanf(ptsp_file, "%s", buf);
            TRACE ( printf("%s ", buf); );
            buf[0]=0;
            fscanf(ptsp_file, "%s", buf);
            TRACE ( printf("%s\n", buf); );
            if ( strcmp("EUC_2D", buf) == 0 )
              {
                distance = round_distance;
              }
            else if ( strcmp("CEIL_2D", buf) == 0 )
              {
                distance = ceil_distance;
              }
            else if ( strcmp("GEO", buf) == 0 )
              {
                distance = geo_distance;
              }
            else if ( strcmp("ATT", buf) == 0 )
              {
                distance = att_distance;
              }
            else
              fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
            strcpy(instance.edge_weight_type, buf);
            buf[0]=0;
          }
        else if( strcmp("EDGE_WEIGHT_TYPE:", buf) == 0 )
          {
            /* set pointer to appropriate distance function; has to be one of
               EUC_2D, CEIL_2D, GEO, or ATT. Everything else fails */
            buf[0]=0;
            fscanf(ptsp_file, "%s", buf);
            TRACE ( printf("%s\n", buf); )
            /*printf("%s\n", buf);
            printf("%s\n", buf);*/
            if ( strcmp("EUC_2D", buf) == 0 )
              {
                distance = round_distance;
              }
            else if ( strcmp("CEIL_2D", buf) == 0 )
              {
                distance = ceil_distance;
              }
            else if ( strcmp("GEO", buf) == 0 )
              {
                distance = geo_distance;
              }
            else if ( strcmp("ATT", buf) == 0 )
              {
                distance = att_distance;
              }
            else
              {
                fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
                exit(1);
              }
            strcpy(instance.edge_weight_type, buf);
            buf[0]=0;
          }
        buf[0]=0;
        fscanf(ptsp_file,"%s", buf);
      }


    if( strcmp("NODE_COORD_SECTION", buf) == 0 )
      {
        TRACE ( printf("found section contaning the node coordinates\n"); )
      }
    else
      {
        fprintf(stderr,"\n\nSome error ocurred finding start of coordinates from tsp file !!\n");
        exit(1);
      }

    if( (nodeptr = malloc(sizeof(struct point) * n)) == NULL )
      exit(EXIT_FAILURE);
    else
      {
        for ( i = 0 ; i < n ; i++ )
          {
            fscanf(ptsp_file,"%ld %lf %lf %lf", &j, &nodeptr[i].x, &nodeptr[i].y, &nodeptr[i].p);
          }
      }

    rewind(ptsp_file);
    fclose(ptsp_file);

    return (nodeptr);

  }



long
int ** compute_distances(void)
/*
      FUNCTION: computes the matrix of all intercity distances
      INPUT:    none
      OUTPUT:   pointer to distance matrix, has to be freed when program stops
*/
{
  long int     i, j;
  long int     **matrix;

  if((matrix = malloc(sizeof(long int) * n * n +
                      sizeof(long int *) * n	 )) == NULL)
    {
      fprintf(stderr,"Out of memory, exit.");
      exit(1);
    }

  for ( i = 0 ; i < n ; i++ )
    {
      matrix[i] = (long int *)(matrix + n) + i*n;
      for ( j = 0  ; j < n ; j++ )
        {
          matrix[i][j] = distance(i,j);
          /*printf("\n the Distance between %d and %d = %d \n",i,j,matrix[i][j]);*/
        }
    }
  return matrix;
}


long int
compute_expected_cost(long int *t)
{
  int no_cities = instance.n,s;
  double expc = 0.0,partial=0.0;
  double *prob_vec=NULL;
  double **Q=NULL;
  int i,r;
  
      prob_vec = (double *) calloc(no_cities,sizeof(double));
      if (!prob_vec)
        error(EXIT_FAILURE,0,"Cannot allocate memory");
      for(i=0; i<no_cities; i++)
        {
          prob_vec[i]= instance.nodeptr[i].p;
        }

      Q = calloc(no_cities,sizeof(double *));
      if(!Q)
        error(EXIT_FAILURE,0,"Cannot allocate memory");
      for(i=0; i < no_cities; i++)
        {
          Q[i] = calloc(no_cities, sizeof(double));
          if(!Q[i])
            error(EXIT_FAILURE,0,"Cannot allocate memory");
        }

      for(i=0; i<no_cities; i++)
        {
          Q[t[i]][t[i]] = 0.0;
          Q[t[i]][t[(i+1)%no_cities]]= 1.0;
		for(s=2; s<no_cities; s++)
                {
                  Q[t[i]][t[(i+s)%no_cities]] = Q[t[i]][t[(i+s-1)%no_cities]]*(1.0 - prob_vec[t[(i+s-1)%no_cities]]);
		  }
        }

      for (i=0; i<n; i++)
        {
          for(r=0; r< n ; r++)
            {
              if(i != r )
                {
                  partial = instance.distance[t[i]][t[r]]
                            * prob_vec[t[i]]
                            * prob_vec[t[r]]
                            * Q[t[i]][t[r]];
                            
                  expc += partial;
                  }
              }
        }

      for(i=0; i < no_cities; i++)
        free(Q[i]);
      free(Q);
      free(prob_vec);

  return (long int)expc;
}

void
printTour( long int *t )
{
  long int   i;

  printf("\n");
  for( i = 0 ; i <= n ; i++ )
    {
      if (!i%25)
        printf("\n");
      printf("%ld ", t[i]);
    }
  printf("\n");

}




const char *argp_program_version = "Evaluator for the PTSP solutions version 1.0";
const char *argp_program_bug_address ="<prasannaprakash@gmail.com>";

static char doc[] =
  "Program to evaluate the solutions of the PTSP.\nThe instance_file sould be in PTSPLIB format.\nThe solution_file should be in a specific format";

static char args_doc[] = "instance_file solution_file";

/* initialise an argp_option struct with the options we except */
static struct argp_option options[] =
    {
      {"time", 't', "TIME",0, "Time up to which the solution cost needs to be computed."
      },
      { 0 }
    };

/* Used by `main' to communicate with `parse_opt'. */
struct arguments
  {
    char *args[2];                /* ARG1 & ARG2 */
    double time;
  };

/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the INPUT argument from `argp_parse', which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = state->input;

  switch (key)
    {
    case 't':
      arguments->time = atof(arg);
      break;
    case ARGP_KEY_ARG:
      if (state->arg_num >= 2)
        /* Too many arguments. */
        argp_usage (state);
      arguments->args[state->arg_num] = arg;
      break;

    case ARGP_KEY_END:
      if (state->arg_num < 2)
        /* Not enough arguments. */
        argp_usage (state);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/* Our argp parser. */
static struct argp argp =
    {
      options, parse_opt, args_doc, doc
    };

int main (int argc, char **argv)
{
  struct arguments arguments;
  char *instance_file_name;
  char *result_file_name;
  char round_str[20];
  int round,j;
  char total_time_str[20];
  double total_time;
  char tour_str[20];
  long int *tour;
  FILE  *result_file;

  /* Default value. */
  arguments.time = 0.0;


  /* Parse arguments; */

  argp_parse (&argp, argc, argv, 0, 0, &arguments);
  instance_file_name=arguments.args[0];
  result_file_name=arguments.args[1];

  instance.nodeptr=read_ptsp(instance_file_name);
  instance.distance=compute_distances();




  if(arguments.time==0.0)
    arguments.time=DBL_MAX;


  tour = (long int *) calloc(instance.n,sizeof(long int));
  if (!tour)
    error(EXIT_FAILURE,0,"Cannot allocate memory");


  result_file = fopen(result_file_name, "r");
  if (result_file == NULL)
    {
      fprintf(stderr,"Can not find result_file, abort\n");
      exit(1);
    }

  while(!feof(result_file))
    {

      fscanf(result_file,"%s %d %s %lf %s",round_str,&round,total_time_str,&total_time,tour_str);


      if(feof(result_file) || (strcmp("Step", round_str)!=0) || (total_time > arguments.time))
        break;

      for( j = 0 ; j < instance.n ; j++ )
        {
          fscanf(result_file,"%ld",&tour[j]);
          if (feof(result_file))
            break;
        }

      if (feof(result_file))
        break;
      printf("%lf\t%ld\n",total_time,compute_expected_cost(tour));

    }


  rewind(result_file);
  fclose(result_file);
  free(tour);
  free(instance.distance);
  free(instance.nodeptr);

  exit(0);
}
