/*
#####################################################################
##########    Estimation-based Local Search for the PTSP    #########
#####################################################################

      Version: 1.0
      File:    heuristics.h
      Authors: Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
      Purpose: Protoypes for heuristic implementations  
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

#include <gsl/gsl_sort_int.h>
#include <gsl/gsl_math.h>

#include "problemdataStructures.h"

#define PI 3.1415

#ifdef __cplusplus
extern "C"
  {
#endif

#ifndef HLS_INCLUDED
# define HLS_INCLUDED

#ifdef LS_SOLUTION_LONGINT
#  define LS_SOLUTION_INT long int
#else
#  define LS_SOLUTION_INT int
#endif
    void printTour( long int *t ) ;
    void compute_tour_length( long int *t ) ;
    void checkTour( long int *t ) ;
    LS_SOLUTION_INT* (*initialSolution)(problem *);
    /* function to compute a solution using the nearest neighbor heuristic*/
    LS_SOLUTION_INT *nearestNeighbor(problem *insPtr);
    /* function to compute a solution using the space filling heuristic*/
    LS_SOLUTION_INT *spaceFilling(problem *insPtr);
    /* function to compute a solution using the radial sort heuristic*/
    LS_SOLUTION_INT *radialSort(problem *insPtr);
    /* function to compute a solution using the nearest insertion heuristic*/
    LS_SOLUTION_INT *nearestInsertion(problem *insPtr);
    /* function to compute a solution using the farthest insertion heuristic*/
    LS_SOLUTION_INT *farthestInsertion(problem *insPtr);
    double radial_angle(double x1, double y1, double x2, double y2);
    int SF_compute_index(int x, int y, long int M, int KMAX);
    /* function to compute quadrent nearest neighbor candidate list*/
    int **computeQuadrantNNLists( problem *insPtr, int nn_ls, int nn_quad);
    void swap(int v[], int v2[], int i, int j);
    void sort(int v[], int v2[], int left, int right);
#endif /* HLS_INCLUDED */


#ifdef __cplusplus

  }
#endif

