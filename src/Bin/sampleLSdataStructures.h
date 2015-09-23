/*
#####################################################################
##########    Estimation-based Local Search for the PTSP    #########
#####################################################################

      Version: 1.0
      File:    sampleLSdataStuctures.h
      Authors: Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
      Purpose: Data structure definition for the local search  
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

#ifndef LS_DATA_STRUCTURES
#define LS_DATA_STRUCTURES


/* Structure describing each city in an apriori solution */
struct LS_city
  {
    int city;
    int dlb;
    double probability;
    double deltaProbability;
    double deltaDashProbability;
    int *realizations;
    int *two_opt_biased_realizations;
    int *two_h_opt_biased_realizations;
    double correction_two_opt[2];
    double correction_two_h_opt;

    //    double *shiftProbabilityVector;
    //    double **correction_vector;
    //    int **geometric_biased_realizations;


    int *neighbors;
    int num_of_zeros;
    int num_of_ones;
    int window_end_node;
    struct LS_city *prev;
    struct LS_city *next;
  };


/* Data structure for describing an apriori solution
   in a local search */
typedef struct
  {
    struct LS_city *array;
    struct LS_city *first;
    int *position_array;
    int no_cities;
    int no_realizations;
    int no_neighbors;
    double value;
    const LS_DISTANCE **distances;
    int generated_realizations;
    int minimum_realizations;
    int maximum_realizations;
    double *delta;
    int *realization_order;
    double alpha;
    int move_status;
    double sum_avg_delta;
    double mean_avg_delta;
    /*importance sampling*/
    int importance_sampling_flag;
    int window_size;
    double no_nodes_inside_window_percentage;

#ifdef LS_EXTRA_STATS_OUTPUT

    int solutions_explored;
    int improvements_made;
    int two_opt_scans_made;
    int two_h_opt_scans_made;
    int samples_used;
#endif

  }
LS_List;

#endif /* LS_DATA_STRUCTURES */
