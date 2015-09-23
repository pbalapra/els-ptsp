/*
#####################################################################
##########    Estimation-based Local Search for the PTSP    #########
#####################################################################

      Version: 1.0
      File:    readFile.c
      Authors: Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
      Purpose: Functions used to read the input file  
      Check:   README and gpl.txt
      Copyright (C) 2008 Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
*/

/*************************************************************************************

    Program's name: els-ptsp

    Estimation-based Local Search for the PTSP 

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
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <values.h>
#include <error.h>
#include "readFile.h"


#define DIST(a,b) insPtr->distance[a][b]
#define TRACE( x )
static double dtrunc (double x)
{
  int k;

  k = (int) x;
  x = (double) k;
  return x;
}

LS_DISTANCE (*distance)(long int i, long int j, problem *insPtr);  /* function pointer */

/*
      FUNCTION: the following four functions implement different ways of 
                computing distances for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
*/

LS_DISTANCE
round_distance (long int i, long int j, problem *insPtr)
/*
      FUNCTION: compute Euclidean distances between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
*/
{
  double xd = insPtr->nodeptr[i].x - insPtr->nodeptr[j].x;
  double yd = insPtr->nodeptr[i].y - insPtr->nodeptr[j].y;
  double r  = sqrt(xd*xd + yd*yd) + 0.5;

  return (LS_DISTANCE) r;

}

LS_DISTANCE
ceil_distance (long int i, long int j, problem *insPtr)
/*
      FUNCTION: compute ceiling distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
*/
{
  double xd = insPtr->nodeptr[i].x - insPtr->nodeptr[j].x;
  double yd = insPtr->nodeptr[i].y - insPtr->nodeptr[j].y;
  double r  = sqrt(xd*xd + yd*yd) + 0.000000001;

  return (LS_DISTANCE)r;
}

LS_DISTANCE
geo_distance (long int i, long int j, problem *insPtr)
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
  double x1 = insPtr->nodeptr[i].x, x2 = insPtr->nodeptr[j].x;
  double y1 = insPtr->nodeptr[i].y, y2 = insPtr->nodeptr[j].y;

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
  return (LS_DISTANCE)dd;

}

LS_DISTANCE
att_distance (long int i, long int j, problem *insPtr)
/*
      FUNCTION: compute ATT distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
*/
{
  double xd = insPtr->nodeptr[i].x - insPtr->nodeptr[j].x;
  double yd = insPtr->nodeptr[i].y - insPtr->nodeptr[j].y;
  /*
  double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
     double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
     */
  double rij = sqrt ((xd * xd + yd * yd) / 10.0);
  double tij = dtrunc (rij);
  long int dij;

  if (tij < rij)
    dij = (int) tij + 1;
  else
    dij = (int) tij;
  return (LS_DISTANCE)dij;
}



struct point
      *read_ptsp(const char *ptsp_file_name, problem *insPtr)
      /*
            FUNCTION: parse and read instance file
            INPUT:    instance name
            OUTPUT:   list of coordinates for all nodes
      	  COMMENTS: Instance files have to be in TSPLIB format, otherwise procedure fails
      */
  {
    FILE         *ptsp_file;
    char         buf[LINE_BUF_LEN];
    long int     i, j,n;
    struct point *nodeptr;


    ptsp_file = fopen(ptsp_file_name, "r");
    if ( ptsp_file == NULL )
      {
        fprintf(stderr,"No instance file specified, abort\n");
        exit(1);
      }
    assert(ptsp_file != NULL);
    /*printf("\nreading tsp-file %s ... \n\n", ptsp_file_name);*/

    fscanf(ptsp_file,"%s", buf);
    while ( strcmp("NODE_COORD_SECTION", buf) != 0 )
      {
        if ( strcmp("NAME", buf) == 0 )
          {
            fscanf(ptsp_file, "%s", buf);
            TRACE ( printf("%s ", buf); )
            fscanf(ptsp_file, "%s", buf);
            /*strcpy(insPtr->name, buf);*/
            TRACE ( printf("%s \n", insPtr->name); )
            buf[0]=0;
          }
        else if ( strcmp("NAME:", buf) == 0 )
          {
            fscanf(ptsp_file, "%s", buf);
            /*strcpy(insPtr->name, buf);*/
            TRACE ( printf("%s \n", insPtr->name); )
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
        else if ( strcmp("TYPE:", buf) == 0 )
          {
            fscanf(ptsp_file, "%s", buf);
            TRACE ( printf("%s\n", buf); )
            if( strcmp("PTSP", buf) != 0 )
              {
                fprintf(stderr,"\n Not a PTSP instance in TSPLIB format !!\n");
                exit(1);
              }
            buf[0]=0;
          }
        else if( strcmp("DIMENSION", buf) == 0 )
          {
            fscanf(ptsp_file, "%s", buf);
            TRACE ( printf("%s ", buf); );
            fscanf(ptsp_file, "%ld", &n);
            insPtr->n = n;

            TRACE ( printf("%ld\n", n); );
            assert ( n > 2 && n < 50000);
            buf[0]=0;
          }
        else if ( strcmp("DIMENSION:", buf) == 0 )
          {
            fscanf(ptsp_file, "%ld", &n);
            insPtr->n = n;
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
            /*strcpy(insPtr->edge_weight_type, buf);*/
            buf[0]=0;
          }
        else if( strcmp("EDGE_WEIGHT_TYPE:", buf) == 0 )
          {
            /* set pointer to appropriate distance function; has to be one of
               EUC_2D, CEIL_2D, GEO, or ATT. Everything else fails */
            buf[0]=0;
            fscanf(ptsp_file, "%s", buf);
            TRACE ( printf("%s\n", buf); )
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
            /*strcpy(insPtr->edge_weight_type, buf);*/
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
            fscanf(ptsp_file,"%ld %lf %lf %lf", &j, &nodeptr[i].x,
                   &nodeptr[i].y,&nodeptr[i].prob);
          }
      }

    fclose(ptsp_file);
    TRACE ( printf("number of cities is %ld\n",n); )
    TRACE ( printf("\n... done\n"); )
    return (nodeptr);

  }


LS_DISTANCE
** compute_distances(problem *insPtr)
/*
      FUNCTION: computes the matrix of all intercity distances
      INPUT:    none
      OUTPUT:   pointer to distance matrix, has to be freed when program stops
*/
{
  long int     i, j,n=insPtr->n;
  LS_DISTANCE     **matrix;


  matrix = malloc(n*sizeof(LS_DISTANCE*));
  if (!matrix)
    error(EXIT_FAILURE,0,"Cannot allocate memory");
  for (i = 0; i < n; i++)
    {
      matrix[i] = malloc(n*sizeof(LS_DISTANCE));
      if (!matrix[i])
        error(EXIT_FAILURE,0,"Cannot allocate memory");
      else
        for ( j = 0  ; j < n ; j++ )
          {
            matrix[i][j] = distance(i,j,insPtr);
          }
    }

  return matrix;
}



