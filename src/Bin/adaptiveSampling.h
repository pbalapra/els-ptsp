/*
#####################################################################
##########    Estimation-based Local Search for the PTSP    #########
#####################################################################

      Version: 1.0
      File:    adaptiveSampling.h
      Authors: Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
      Purpose: Function prototypes for implementing adaptive sampling and improtance sampling  
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

double
(*delta_evaluation)(LS_List *solPtr,
                    int edge0fst, int edge0snd, int node,
                    int edge1fst, int edge1snd, gsl_rng *r, int opt);  /* function pointer */

double
LSA_2opt_delta_sample_estimate(LS_List *solPtr,
                               int edge0fst, int edge0snd,
                               int edge1fst, int edge1snd, int realization_index,int importance_sampling);
double
LSA_2hopt_delta_sample_estimate(LS_List *solPtr, int edge0fst, int edge0snd, int node,
                                int edge1fst, int edge1snd, int realization_index,int importance_sampling);

void
LSA_solution_add_realization(LS_List solution,
                             int realization_number, gsl_rng *r);

double
LS_delta(LS_List *solPtr,
         int edge0fst, int edge0snd, int node,
         int edge1fst, int edge1snd,
         gsl_rng *r, int opt);

double
LSA_delta_adaptive_sample(LS_List *solPtr,
                          int edge0fst, int edge0snd, int node,
                          int edge1fst, int edge1snd, gsl_rng *r, int opt);

double
LSA_2opt_delta_sample_estimate_window(LS_List *solPtr,
                                      int edge0fst, int edge0snd,
                                      int edge1fst, int edge1snd, int realization_index,int importance_sampling);
double
LSA_2hopt_delta_sample_estimate_window(LS_List *solPtr, int edge0fst, int edge0snd, int node,
                                       int edge1fst, int edge1snd, int realization_index,int importance_sampling);
