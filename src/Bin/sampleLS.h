/*
#####################################################################
##########    Estimation-based Local Search for the PTSP    #########
#####################################################################

      Version: 1.0
      File:    sampleLs.h
      Authors: Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
      Purpose: Function prototypes for neighborhood exploration  
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

#include <gsl/gsl_rng.h>
#include "problemdataStructures.h"

#ifdef __cplusplus
extern "C"
  {
#endif


#ifndef LS_INCLUDED
# define LS_INCLUDED

#ifdef LS_SOLUTION_LONGINT
#  define LS_SOLUTION_INT long int
#else
#  define LS_SOLUTION_INT int
#endif

#if defined(LS_DISTANCE_LONGINT) && defined(LS_DISTANCE_DOUBLE)
#error "-DLS_DISTANCE_LONGINT conflicts with -DLS_DISTANCE_DOUBLE"
#endif

#ifdef LS_DISTANCE_LONGINT
#  define LS_DISTANCE long int
#elif LS_DISTANCE_DOUBLE
#  define LS_DISTANCE double
#else
#  define LS_DISTANCE int
#endif

    /* include definition of data structures for local search */
#include "sampleLSdataStructures.h"

    /* include definition of some auxiliary functions */
#include "sampleLSauxiliaryFunctions.h"

    /* Allocate a data structure for a local search. The parameter D is
       the pointer to the matrix containing the distances between
       cities. The matrix D IS_NOT modified: To insure this, it is
       immediatelly casted to `const LS_DISTANCE'.  This function is
       declared here as `extern inline'. A static version of the function
       is given in sampleLScommon.c for the cases in which the compiler
       does not inline. */
    extern inline LS_List
      LS_solution_allocate(int no_cities, int no_realizations,
                           const double *prob_vec, LS_DISTANCE **D,
                           double alpha, int importance_sampling,
                           float deltaProb, float deltaDashProb,
                           float window_size_percent,float nodes_percent)
      {
        return LS_solution_allocate_aux(no_cities,no_realizations,prob_vec,(const LS_DISTANCE **) D,
                                        alpha,importance_sampling,
                                        (double) deltaProb, (double) deltaDashProb,window_size_percent,nodes_percent);
      }


    /* In the following, we assume that apriori_solution is a vector of
       integer (or long integers if LS_SOLUTION_LONGINT is defined)
       describing an apriori solution. */
    void
    LS_solution_set(LS_List *solPtr, const LS_SOLUTION_INT *apriori_solution);

    void
    LS_solution_set_with_computed_value(LS_List *solPtr,
                                        const LS_SOLUTION_INT *apriori_solution);
    void
    LS_solution_set_with_given_value(LS_List *solPtr, double value,
                                     const LS_SOLUTION_INT *apriori_solution);

    /* The following function is to be called no_realizations times. */
    void
    LS_solution_add_realization(LS_List solution,
                                int realization_number,
                                const double *ran_num);

    /* The following functions adds a number of realizations equal to
       "solution.no_realizations". The probabilities used for generating
       the realizations are those store in the variable "solution". */
    void
    LS_resample_realizations(LS_List solution, gsl_rng *r, int sampling_type);



    /* A single step of a 2-exchange first improvement local search with
      nearest neighbor lists and don't look bits */
    double
    LS_2nndlbfls_step(LS_List *solPtr, const int *order, gsl_rng *r, int sampling_type);

    /* Iterate 2nndlbfls */
    void
    LS_2nndlbfls(LS_List *solPtr, gsl_rng *r, double time, int verbose, int sampling_type);

    /* Iterate 2nndlbfls K times */
    void
    LS_2nndlbfls_times(int K,LS_List *solPtr, gsl_rng *r,
                       double time, int verbose,int sampling_type);


    /* A single step of a 2.5-exchange first improvement local search with
       nearest neighbor lists and don't look bits */
    double
    LS_2hnndlbfls_step(LS_List *solPtr, const int *order, gsl_rng *r, int sampling_type);

    /* Iterate 2hnndlbfls */
    void
    LS_2hnndlbfls(LS_List *solPtr, gsl_rng *r, double time, int verbose, int sampling_type);

    /* Iterate 2hnndlbfls K times */
    void
    LS_2hnndlbfls_times(int K,LS_List *solPtr, gsl_rng *r,
                        double time, int verbose, int sampling_type);


    /* Iterate 2nndlbfls, resample realizations at each step. */
    void
    LS_2nndlbfls_resample(LS_List *solPtr,
                          gsl_rng *r, double time, int verbose, int sampling_type);

    /* Iterate 2nndlbfls K times, resample realizations at each step. */
    void
    LS_2nndlbfls_times_resample(int K, LS_List *solPtr,
                                gsl_rng *r, double time, int verbose,int sampling_type);

    /* Iterate 2hnndlbfls, resample realizations at each step. */
    void
    LS_2hnndlbfls_resample(LS_List *solPtr,
                           gsl_rng *r, double time, int verbose, int sampling_type);

    /* Iterate 2hnndlbfls K times, resample realizations at each step. */
    void
    LS_2hnndlbfls_times_resample(int K,LS_List *solPtr,
                                 gsl_rng *r, double time, int verbose,int sampling_type);

    /* Reset don't look bits of a solution */
    void
    LS_reset_dlb(LS_List *solPtr);

    /* Compute the value of the apriori solution */
    double
    LS_solution_return_value(LS_List solution);

    /* Compute and set the value of the apriori solution */
    double
    LS_solution_compute_and_set_value(LS_List *solPtr);

    /* Writes the currrent solution into the vector pointed by
       apriori_solution. Returns the value of the current solution. */
    double
    LS_solution_get(LS_List solution, LS_SOLUTION_INT *apriori_solution);

    /* Dispose of a solution list */
    void
    LS_solution_free(LS_List *solPtr);

    /* Print a solution... just for debugging purposes */
    void
    LS_solution_print(LS_List solution);
    void
    LS_solution_print_all(LS_List solution);

    /* Log results to stdout */
    void
    LS_solution_log(LS_List solution, int i);


    /* Sort neighbors of each city. The parameter `nn' is the number of
       neighbors. IMPORTANT: a city is the nearest neighbor of itself!!!
       So it always appears in position 0 in its nearest neighbors
       list. */
    void
    LS_solution_sort_neighbors(LS_List *solPtr, int nn);

    void
    LS_solution_sort_quad_neighbors(problem *insPtr, LS_List *solPtr,int nn);

    int*
    LS_allocate_sort_quad_neighbors(problem *insPtr,const LS_DISTANCE **d,int node,int no_cities,int nn);

    void
    LS_swap(int v[], int v2[], int i, int j);

    void
    LS_sort(int v[], int v2[], int left, int right);

    double
    LS_compute_expected_cost(int no_cities,long int* apriori_solution,
                             LS_DISTANCE **D, double *prob_vec,int homoflag);

    int
    LS_move_check_position_index(LS_List *solPtr);

#endif /* LS_INCLUDED */


#ifdef __cplusplus

  }
#endif
