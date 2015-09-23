/*
#####################################################################
##########    Estimation-based Local Search for the PTSP    #########
#####################################################################

      Version: 1.0
      File:    heuristics.c
      Authors: Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
      Purpose: Functions for implementing several heuristics
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
#include <error.h>
#include <values.h>

#include "heuristics.h"


#  include <gsl/gsl_sort_double.h>

#define TRUE 1
#define FALSE 0
#define INFTY    LONG_MAX

#ifdef ARRAYTEST
#define DIST(a,b) DISTa[a][b]
#define BIG 0   /* book is wrong */
static int DISTa[6][6] = {
                           {BIG,   3,      93,     13,     33,     9},
                           {4,     BIG,    77,     42,     21,     16},
                           {45,    17,     BIG,    36,     16,     28},
                           {39,    90,     80,     BIG,    56,     7},
                           {28,    46,     88,     33,     BIG,    25},
                           {3,     88,     18,     46,     92,     BIG}};
#else
/* warning: this macro cannot be called twice in the same arithmetic
* expression*/
#define DIST(a,b) insPtr->distance[a][b]
#endif

void printTour( long int *t )
/*
      FUNCTION:       print the tour *t
      INPUT:          pointer to a tour
      OUTPUT:         none
*/
{
  /*   long int   i;

     printf("\n");
     for( i = 0 ; i <= n ; i++ ) {
     if (!i%25) printf("\n");
     printf("%ld ", t[i]);
      }
     printf("\n");
  */
}



void compute_tour_length( long int *t )
/*
      FUNCTION: compute the tour length of tour t
      INPUT:    pointer to tour t
      OUTPUT:   tour length of tour t
*/
{
  /*   int      i;
     long int tour_length = 0;

  for ( i = 0 ; i < n ; i++ ) {
  tour_length += instance.distance[t[i]][t[i+1]];
     }


     printf("The tour lenght =%ld\n",tour_length);
  */

}

void checkTour( long int *t )
/*
      FUNCTION:       make a simple check whether tour *t can be feasible
      INPUT:          pointer to a tour
      OUTPUT:         none
*/
{
  /* long int   i, sum=0;

   for( i = 0 ; i < n ; i++ ) {
  sum += t[i];
   }
   if ( sum != (n-1) * n / 2 ) {
  fprintf(stderr,"Next tour must be flawed !!\n");
  printTour( t );
  exit(1);
   }*/
}

LS_SOLUTION_INT *nearestNeighbor(problem *insPtr)
{
  int i, j,n=insPtr->n;
  LS_SOLUTION_INT  *tour;
  unsigned long int tourlength=0;
  int *visited;
  long int this_city;
  long int closest_city=-1;
  double closest_distance;
  int start_city=-1;

  tour=malloc(n*sizeof(LS_SOLUTION_INT));
  if (!tour)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  visited=malloc(n*sizeof(int));
  if (!visited)
    error(EXIT_FAILURE,0,"Cannot allocate memory");


  /* initialize unvisited cities */
  for (i = 0; i < n; i++)
    {
      visited[i] = FALSE;
    }

  /* choose number_of_cities as starting point */
  this_city = 0;
  tour[0]=0;
  visited[this_city] = TRUE;
  if(start_city==-1)
    start_city=this_city;
  /* main loop of nearest neighbor heuristic */
  for (i = 1; i < n; i++)
    {

      /* find nearest unvisited city to this city */
      closest_distance = INFTY;
      for (j = 0; j < n; j++)
        {
          if (!visited[j])
            {
              /*printf("From %ld to %d the tour lenght =%ld\n",this_city,j,insPtr.distance[this_city][j]);*/
              if (insPtr->distance[this_city][j] < closest_distance)
                {
                  closest_distance = insPtr->distance[this_city][j];
                  closest_city = j;
                }
            }
        }
      /* report closest city */





      tourlength+=insPtr->distance[this_city][closest_city];
      tour[i]=closest_city;
      visited[closest_city] = TRUE;
      this_city = closest_city;


    }

  /* finish tour by returning to start */

  tourlength+=insPtr->distance[this_city][start_city];
  //    tour[n]=0;
  checkTour(tour);
  /*printf("\nNearest Neighbour\n");
  printf("Tourlength %5ld \t",tourlength);*/
  /*printf("TOUR\t");
  printTour(tour);
  printf("\n");*/
  compute_tour_length(tour);
  free(visited);
  return tour;
  /*free(tour);*/
}



LS_SOLUTION_INT *nearestInsertion(problem *insPtr)


{
  int end1=-1, end2=-1,n=insPtr->n;
  int s=0, nih=1;
  int i, j, argmaxcost=-1;
  int index, nextindex;
  long int inscost=0, newcost=0, total_cost=0, maxcost=0;
  int *cycle; /* cycle[i]=next node after node i in tour; -1 if not yet in */
  long int *dist;  /* dist[i] =dist from any node in tour to node i not in tour */
  long int  *route;

  assert(n > 0);

  route=malloc(n*sizeof(LS_SOLUTION_INT));
  if (!route)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  cycle = malloc(n*sizeof(long int));
  if (!cycle)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  dist  = malloc(n*sizeof(long int));
  if (!dist)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  /* No nodes in tour yet */
  for (i=0; i<n; i++)
    cycle[i] = -1;
  /* Start with a one-node tour starting and ending at node s */
  cycle[s] = s;
  total_cost = 0;
  /* Find distance of all nodes from closest node in tour */
  for (i=0; i<n; i++)
    {
      dist[i] = DIST(s,i);
    }

  /* Until all nodes have been added to the tour */
  for (i=0; i<n-1; i++)
    {

      /* Find node which is farthest from (nearest to) the tour */
      maxcost = nih ? MAXINT : -MAXINT;
      argmaxcost = -1;
      for (j=0; j<n; j++)
        {
          if (cycle[j] < 0)
            {
              if ((!nih && (dist[j] > maxcost)) || (nih&&(dist[j]<maxcost)))
                {
                  maxcost = dist[j];
                  argmaxcost = j;
                }
            }
        }

      assert(argmaxcost >= 0 && argmaxcost < n);

      /* Farthest node is argmaxcost, distance is maxcost. */
      /* Now find cheapest place to insert it into the tour. */
      inscost = MAXINT;
      index = s;
      for (j=0; j<=i; j++)
        {
          nextindex = cycle[index];
          assert(nextindex >= 0 && nextindex < n);
          newcost = DIST(index,argmaxcost);
          newcost += DIST(argmaxcost, nextindex);
          newcost -= DIST(index, nextindex);
          if (newcost < inscost)
            {
              inscost = newcost;
              end1 = index;
              end2 = nextindex;
            }
          index = nextindex;
        }

      /* Cheapest place is after end1 and before end2.  Insert it there. */
      assert(cycle[end1] == end2);
      cycle[end1] = argmaxcost;
      cycle[argmaxcost] = end2;
      total_cost += inscost;


      /* Now that there's a new node in the tour, update the array of
       * distances from other nodes to the tour.
       */
      for (j=0; j<n; j++)
        {
          /* If node j not in tour, */
          if (cycle[j] < 0)
            {
              /* If closer to the new node than to any other node in the tour,
               * record the new distance.
               */
              newcost = DIST(argmaxcost, j);
              if (newcost < dist[j])
                dist[j] = newcost;
            }
        }


    }

  /* Find longest move in tour. */
  index = s;
  maxcost = -MAXINT;
  for (i=0; i<n; i++)
    {
      newcost = DIST(index, cycle[index]);
      if (newcost > maxcost)
        {
          maxcost = newcost;
          argmaxcost = index;
        }
      index = cycle[index];
    }


  /* Return the tour to the caller as an array of nodes to visit
   * rather than a map {current node, next node} as it is in dist[].
   * Start tour such that longest move in tour is (route[n-1],route[0]).
   */
#ifdef ARRAYTEST

  maxcost = 0;
  index = s;
#else

  index = cycle[argmaxcost];
#endif

  for (i=0; i<n; i++)
    {
      route[i] = index;
      index = cycle[index];
    }
  //    route[n]=route[0];

//  printf("Nearest inseriton\n");

  checkTour(route);
  /*printf("Tourlength = %5ld \t",total_cost - maxcost+DIST(route[n-1],route[0]));*/
  /*printf("TOUR\t");
  printTour(route);
  printf("\n");*/
  compute_tour_length(route);

  free(cycle);
  free(dist);

  return route;
  /* Return length of one-way trip */

}



LS_SOLUTION_INT *farthestInsertion(problem *insPtr)


{
  int end1=-1, end2=-1,n=insPtr->n;
  int s=0, nih=0;
  int i, j, argmaxcost=-1;
  int index, nextindex;
  long int inscost=0, newcost=0, total_cost=0, maxcost=0;
  int *cycle; /* cycle[i]=next node after node i in tour; -1 if not yet in */
  long int *dist;  /* dist[i] =dist from any node in tour to node i not in tour */
  long int  *route;

  assert(n > 0);


  route=malloc(n*sizeof(LS_SOLUTION_INT));
  if (!route)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  cycle = malloc(n*sizeof(long int));
  if (!cycle)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  dist  = malloc(n*sizeof(long int));
  if (!dist)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  /* No nodes in tour yet */
  for (i=0; i<n; i++)
    cycle[i] = -1;
  /* Start with a one-node tour starting and ending at node s */
  cycle[s] = s;
  total_cost = 0;
  /* Find distance of all nodes from closest node in tour */
  for (i=0; i<n; i++)
    {
      dist[i] = DIST(s,i);

    }

  /* Until all nodes have been added to the tour */
  for (i=0; i<n-1; i++)
    {

      /* Find node which is farthest from (nearest to) the tour */
      maxcost = nih ? MAXINT : -MAXINT;
      argmaxcost = -1;
      for (j=0; j<n; j++)
        {
          if (cycle[j] < 0)
            {
              if ((!nih && (dist[j] > maxcost)) || (nih&&(dist[j]<maxcost)))
                {
                  maxcost = dist[j];
                  argmaxcost = j;
                }
            }
        }



      assert(argmaxcost >= 0 && argmaxcost < n);

      /* Farthest node is argmaxcost, distance is maxcost. */
      /* Now find cheapest place to insert it into the tour. */
      inscost = MAXINT;
      index = s;
      for (j=0; j<=i; j++)
        {
          nextindex = cycle[index];
          assert(nextindex >= 0 && nextindex < n);
          newcost = DIST(index,argmaxcost);
          newcost += DIST(argmaxcost, nextindex);
          newcost -= DIST(index, nextindex);
          if (newcost < inscost)
            {
              inscost = newcost;
              end1 = index;
              end2 = nextindex;
            }
          index = nextindex;
        }

      /* Cheapest place is after end1 and before end2.  Insert it there. */
      assert(cycle[end1] == end2);
      cycle[end1] = argmaxcost;
      cycle[argmaxcost] = end2;
      total_cost += inscost;


      /* Now that there's a new node in the tour, update the array of
       * distances from other nodes to the tour.
       */
      for (j=0; j<n; j++)
        {
          /* If node j not in tour, */
          if (cycle[j] < 0)
            {
              /* If closer to the new node than to any other node in the tour,
               * record the new distance.
               */
              newcost = DIST(argmaxcost, j);
              if (newcost < dist[j])
                dist[j] = newcost;
            }
        }


    }

  /* Find longest move in tour. */
  index = s;
  maxcost = -MAXINT;
  for (i=0; i<n; i++)
    {
      newcost = DIST(index, cycle[index]);
      if (newcost > maxcost)
        {
          maxcost = newcost;
          argmaxcost = index;
        }
      index = cycle[index];
    }


  /* Return the tour to the caller as an array of nodes to visit
   * rather than a map {current node, next node} as it is in dist[].
   * Start tour such that longest move in tour is (route[n-1],route[0]).
   */
#ifdef ARRAYTEST

  maxcost = 0;
  index = s;
#else

  index = cycle[argmaxcost];
#endif

  for (i=0; i<n; i++)
    {
      route[i] = index;
      index = cycle[index];
    }
  //    route[n]=route[0];

/*
  if(!nih)
    printf("Fartherst inseriton\n");
*/    
  checkTour(route);
  /*printf("Tourlength = %5ld \t",total_cost - maxcost+DIST(route[n-1],route[0]));*/
  /*printf("TOUR\t");
  printTour(route);
  printf("\n");*/
  compute_tour_length(route);

  free(cycle);
  free(dist);

  return route;

}




double radial_angle(double x1, double y1, double x2, double y2)
{
  float dx = x2-x1;
  float dy = y2-y1;
  double angle=0.0;

  // Calculate angle
  if (dx == 0.0)
    {
      if (dy == 0.0)
        angle = 0.0;
      else if (dy > 0.0)
        angle = PI / 2.0;
      else
        angle = PI * 3.0 / 2.0;
    }
  else if (dy == 0.0)
    {
      if  (dx > 0.0)
        angle = 0.0;
      else
        angle = PI;
    }
  else
    {
      if  (dx < 0.0)
        angle = atan(dy/dx) + PI;
      else if (dy < 0.0)
        angle = atan(dy/dx) + (2*PI);
      else
        angle = atan(dy/dx);
    }

  // Convert to degrees
  angle = angle * 180 / PI;

  // Return
  return angle;
}




int SF_compute_index(int x, int y, long int M, int KMAX)
{
  int index=0,K=1,temp;

  if(x>y)
    {
      index=1;
      x=M-x;
      y=M-y;
    }

  while(K<KMAX)
    {
      index=index+1;
      K=K+1;
      if( (x+y) > M)
        {
          index=index+1;
          temp=M-y;
          y=x;
          x=temp;
        }

      if(K<KMAX)
        {
          index=index+1;
          K=K+1;
          x=x+x;
          y=y+y;
          if(y>M)
            {
              index=index+1;
              temp=y-M;
              y=M-x;
              x=temp;
            }
        }
    }

  return index;
}

LS_SOLUTION_INT *spaceFilling(problem *insPtr)
{

  int *sierpinski_index;
  long int* sorted_tour;
  int i,n=insPtr->n;

  sierpinski_index= malloc(n*sizeof(double));
  if (!sierpinski_index)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  sorted_tour = malloc(n*sizeof(LS_SOLUTION_INT));
  if (!sorted_tour)
    error(EXIT_FAILURE,0,"Cannot allocate memory");


  for(i=0;i<n;i++)
    {
      sierpinski_index[i]=SF_compute_index(floor(insPtr->nodeptr[i].x),
                                           floor(insPtr->nodeptr[i].y), 1,n);
      //printf("index of the city %d = %d\n",i,sierpinski_index[i]);
    }


  gsl_sort_int_smallest_index((size_t*)sorted_tour,n,sierpinski_index,1,n);

  printTour(sorted_tour) ;

  free(sierpinski_index);
  return sorted_tour;

}


LS_SOLUTION_INT *radialSort(problem *insPtr)
{

  double sum_of_x=0.0,sum_of_y=0.0,center_of_x=0.0,center_of_y=0.0;
  double *angle;
  long int* sorted_tour;
  int i,n=insPtr->n;

  angle = malloc(n*sizeof(double));
  if (!angle)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  sorted_tour = malloc(n*sizeof(LS_SOLUTION_INT));
  if (!sorted_tour)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  for(i=0;i<n;i++)
    {
      sum_of_x+=insPtr->nodeptr[i].x;
      sum_of_y+=insPtr->nodeptr[i].y;
    }

  center_of_x=sum_of_x/(double)n;
  center_of_y=sum_of_y/(double)n;


  for(i=0;i<n;i++)
    {
      angle[i]=radial_angle(center_of_x, center_of_y, insPtr->nodeptr[i].x,
                            insPtr->nodeptr[i].y);
      //printf("angle of city %d = %f\n",i,angle[i]);
    }

  gsl_sort_smallest_index((size_t*)sorted_tour,n,angle,1,n);

  printTour(sorted_tour) ;

  free(angle);

  return sorted_tour;

}
