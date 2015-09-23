/*
#####################################################################
##########    Estimation-based Local Search for the PTSP    #########
#####################################################################

      Version: 1.0
      File:    adaptiveSampling.c
      Authors: Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
      Purpose: Functions for implementing adaptive sampling and improtance sampling  
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
#include <string.h>
#include <stdio.h>
#include <error.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_cdf.h>
#include <search.h>

#ifdef LS_DISTANCE_LONGINT
#  include <gsl/gsl_sort_long.h>
#elif LS_DISTANCE_DOUBLE
#  include <gsl/gsl_sort_double.h>
#else
#  include <gsl/gsl_sort_int.h>
#endif



#include "sampleLS.h"
#include "adaptiveSampling.h"
#include "stopwatch.h"
#define TRUE	1
#define FALSE	0

#include "statsTables.h"


inline void
LSA_initialize_zero_one_window(LS_List *solPtr,int node)
{
  solPtr->array[node].num_of_ones=0;
  solPtr->array[node].num_of_zeros=0;
  solPtr->array[node].window_end_node=-1;
}

inline int
LSA_segment_length(LS_List *solPtr,int node1, int node2)
{
  int n=solPtr->no_cities,k=0,i=-1,j=-1;
  i=solPtr->position_array[node1];
  j=solPtr->position_array[node2];
  //printf("pos:%d\t%d\n",i,j);

  if(j==i)
    k=0;
  else if (j>i)
    k = j-i;
  else
    k = j+n-i;
  return k;
}


double
LSA_2opt_delta_sample_estimate(LS_List *solPtr,
                               int edge0fst, int edge0snd,
                               int edge1fst, int edge1snd,int index_realization, int importance_sampling)
{
  int cover0fst, cover0snd, cover1fst, cover1snd;
  long int delta = 0;
  int j;
  struct LS_city *this;
  const LS_DISTANCE **D = solPtr->distances;
  int window_length=0;
  int imp_flag;
  double correction=1.0;

  for (j=index_realization; j<index_realization+1; j++)
    {
      cover0fst = edge0fst;
      cover0snd = edge0snd;

      /* Search to the right of insertion point 0 */
      this = &solPtr->array[cover0snd];


      imp_flag=1;

      while (this->city != edge1snd)
        {
          if(imp_flag && importance_sampling)
            {
              if (this->two_opt_biased_realizations[j])
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[1];
                  break;
                }
              else
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[0];
                  this = this->next;
                }
            }
          else
            {
              imp_flag=0;
              if (this->realizations[j])
                break;
              this = this->next;
            }
#ifdef LS_EXTRA_STATS_OUTPUT
          solPtr->two_opt_scans_made++;
#endif

        }


      if (this->city == edge1snd)
        return 0.0;

      cover0snd = this->city;

      /* Search to the left of insertion point 0 */
      this = &solPtr->array[cover0fst];
      window_length=LSA_segment_length(solPtr,edge1fst,this->city);

      imp_flag=1;

      while (this->city != edge1fst)
        {
          if(imp_flag && importance_sampling )
            {
              if (this->two_opt_biased_realizations[j])
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[1];
                  break;
                }
              else
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[0];
                  this = this->prev;
                }
            }
          else
            {
              imp_flag=0;
              if (this->realizations[j])
                break;
              this = this->prev;
            }
#ifdef LS_EXTRA_STATS_OUTPUT
          solPtr->two_opt_scans_made++;
#endif

        }

      if (this->city == edge1fst)
        return 0.0;

      cover0fst = this->city;

      cover1fst = edge1fst;
      cover1snd = edge1snd;

      /* Search to the right of insertion point 1 */
      this = &solPtr->array[cover1snd];

      imp_flag=1;

      while (this->city != edge0snd)
        {
          if(imp_flag && importance_sampling  )
            {
              if (this->two_opt_biased_realizations[j])
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[1];
                  break;
                }
              else
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[0];
                  this = this->next;
                }
            }
          else
            {
              imp_flag=0;
              if (this->realizations[j])
                break;
              this = this->next;
            }
#ifdef LS_EXTRA_STATS_OUTPUT
          solPtr->two_opt_scans_made++;
#endif

        }


      if (this->city == edge0snd)
        return 0.0;
      cover1snd = this->city;

      /* Search to the left of insertion point 1 */
      this = &solPtr->array[cover1fst];

      imp_flag=1;

      while (this->city != edge0fst)
        {

          if( imp_flag && importance_sampling )
            {
              if (this->two_opt_biased_realizations[j])
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[1];
                  break;
                }
              else
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[0];
                  this = this->prev;
                }

            }
          else
            {
              imp_flag=0;
              if (this->realizations[j])
                break;
              this = this->prev;
            }

#ifdef LS_EXTRA_STATS_OUTPUT
          solPtr->two_opt_scans_made++;
#endif

        }

      if (this->city == edge0fst)
        return 0.0;
      cover1fst = this->city;

      delta =-D[cover0fst][cover0snd]
             -D[cover1fst][cover1snd]
             +D[cover0fst][cover1fst]
             +D[cover1snd][cover0snd];
    }


  return ((double)delta*correction);


}

double
LSA_2hopt_delta_sample_estimate(LS_List *solPtr,
                                int edge0fst, int edge0snd, int node,
                                int edge1fst, int edge1snd, int index_realization,int importance_sampling)
{
  int cover0fst, cover0snd, cover1fst, cover1snd;
  long int sum_delta = 0;
  int j;
  struct LS_city *this=NULL;
  const LS_DISTANCE **D = solPtr->distances;

  double correction=1.0;
  int imp_flag;


  for (j=index_realization; j<index_realization+1; j++)
    {
      if(importance_sampling)
        {
          if (!(solPtr->array[node].two_h_opt_biased_realizations[j]))
            {
              return 0.0;
            }

          correction*=solPtr->array[node].correction_two_h_opt;
        }
      else
        {
          if (!solPtr->array[node].realizations[j])
            {
              return 0.0;
            }
        }


      cover0fst = edge0fst;
      cover0snd = edge0snd;


      if(importance_sampling==4)
        importance_sampling=0;

      /* Search to the right of insertion point 0 */
      this = &solPtr->array[cover0snd];

      imp_flag=1;

      while (this->city != edge1snd)
        {
          if(imp_flag && importance_sampling)
            {
              if (this->two_opt_biased_realizations[j])
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[1];
                  break;
                }
              else
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[0];
                  this = this->next;
                }
            }
          else
            {
              imp_flag=0;
              if (this->realizations[j])
                break;
              this = this->next;
            }
#ifdef LS_EXTRA_STATS_OUTPUT
          solPtr->two_opt_scans_made++;
#endif

        }




      if (this->city == edge1snd)
        return 0.0;

      cover0snd = this->city;



      /* Search to the left of insertion point 0 */
      this = &solPtr->array[cover0fst];
      imp_flag=1;

      while (this->city != edge1fst)
        {
          if(imp_flag && importance_sampling)
            {
              if (this->two_opt_biased_realizations[j])
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[1];
                  break;
                }
              else
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[0];
                  this = this->prev;
                }
            }
          else
            {
              imp_flag=0;
              if (this->realizations[j])
                break;
              this = this->prev;
            }
#ifdef LS_EXTRA_STATS_OUTPUT
          solPtr->two_opt_scans_made++;
#endif

        }


      if (this->city == edge1fst)
        return 0.0;

      cover0fst = this->city;

      cover1fst = edge1fst;
      cover1snd = edge1snd;


      /* Search to the right of insertion point 1 */
      this = &solPtr->array[cover1snd];

      imp_flag=1;

      while (this->city != edge0snd)
        {
          if(imp_flag && importance_sampling )
            {
              if (this->two_opt_biased_realizations[j])
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[1];
                  break;
                }
              else
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[0];
                  this = this->next;
                }
            }
          else
            {
              imp_flag=0;
              if (this->realizations[j])
                break;
              this = this->next;
            }
#ifdef LS_EXTRA_STATS_OUTPUT
          solPtr->two_opt_scans_made++;
#endif

        }




      if (this->city == edge0snd)
        return 0.0;
      cover1snd = this->city;



      /* Search to the left of insertion point 1 */
      this = &solPtr->array[cover1fst];

      imp_flag=1;

      while (this->city != edge0fst)
        {

          if( imp_flag && importance_sampling )
            {
              if (this->two_opt_biased_realizations[j])
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[1];
                  break;
                }
              else
                {
                  imp_flag=0;
                  correction*=this->correction_two_opt[0];
                  this = this->prev;
                }

            }
          else
            {
              imp_flag=0;
              if (this->realizations[j])
                break;
              this = this->prev;
            }

#ifdef LS_EXTRA_STATS_OUTPUT
          solPtr->two_opt_scans_made++;
#endif

        }

      if (this->city == edge0fst)
        return 0.0;
      cover1fst = this->city;


      sum_delta +=-D[cover0fst][cover0snd]
                  -D[cover1fst][node]
                  -D[node][cover1snd]
                  +D[cover0fst][node]
                  +D[node][cover0snd]
                  +D[cover1fst][cover1snd];
    }


  return ((double)sum_delta*correction);

}

void
LSA_solution_add_realization(LS_List solution,
                             int realization_number,
                             gsl_rng *r)
{
  int i,j;
  int *realization;

  realization = malloc(solution.no_cities*sizeof(int));
  for (j=0; j<solution.no_cities; j++)
    realization[j] = (gsl_rng_uniform(r)<=solution.array[j].probability);

  for (i=0;i<solution.no_cities;i++)
    {
      solution.array[i].realizations[realization_number] = realization[i];
    }
  free(realization);
}


double
LS_delta(LS_List *solPtr,
         int edge0fst, int edge0snd, int node,
         int edge1fst, int edge1snd,
         gsl_rng *r, int opt)
{
  double sum_delta = 0.0,delta=0.0;
  int index_realization=0;
  solPtr->move_status=FALSE;

  //printf("solPtr->no_realizations=%d\n",solPtr->no_realizations);
  //printf("solPtr->importance_sampling_flag=%d\n",solPtr->importance_sampling_flag);
  while(index_realization < solPtr->no_realizations)
    {
#ifdef LS_EXTRA_STATS_OUTPUT
      solPtr->samples_used++;
#endif

      if(opt==0)
        delta= LSA_2opt_delta_sample_estimate(solPtr,edge0fst, edge0snd,
                                              edge1fst, edge1snd,index_realization, 0);
      else
        delta= LSA_2hopt_delta_sample_estimate(solPtr,edge0fst, edge0snd,node,
                                               edge1fst, edge1snd,index_realization, 0);
      sum_delta +=delta;
      index_realization++;
    }


  if (sum_delta < 0.0)
    {
      solPtr->move_status=TRUE;
    }

#ifdef LS_EXTRA_STATS_OUTPUT
  solPtr->solutions_explored++;
#endif

  return (sum_delta/(double)solPtr->no_realizations);
}


double
LSA_delta_adaptive_sample(LS_List *solPtr,
                          int edge0fst, int edge0snd, int node,int edge1fst, int edge1snd,
                          gsl_rng *r, int opt)
{
  long int sum_delta = 0;
  double delta=0.0;
  double mean_delta=0.0;
  double sd_delta=0.0;
  double sd_estimator_delta=0.0;
  double sum_delta_sqr=0.0;
  int index_realization=0;
  int realization=0;
  double chebyshev_k=0.0;
  int alpha_index=0;
  double alpha=solPtr->alpha;

  if(!gsl_fcmp(alpha,0.10,0.0000001))
    {
      alpha_index=0;
    }
  else if(!gsl_fcmp(alpha,0.05,0.0000001))
    {
      alpha_index=1;
    }
  else if(!gsl_fcmp(alpha,0.02,0.0000001))
    {
      alpha_index=2;
    }
  else if(!gsl_fcmp(alpha,0.01,0.0000001))
    {
      alpha_index=3;
    }
  else
    {
      printf("wrong alpha type\n");
      exit(0);
    }

  solPtr->move_status=FALSE;


  for(index_realization=0; index_realization< solPtr->maximum_realizations;index_realization++)
    {

#ifdef LS_EXTRA_STATS_OUTPUT
      solPtr->samples_used++;
#endif

      realization=solPtr->realization_order[index_realization];
      if(opt==0)
        {
          if(solPtr->importance_sampling_flag==0)
            delta=LSA_2opt_delta_sample_estimate(solPtr,edge0fst, edge0snd,
                                                 edge1fst, edge1snd, realization,0);
          else
            delta=LSA_2opt_delta_sample_estimate_window(solPtr,edge0fst, edge0snd,
                  edge1fst, edge1snd, realization,solPtr->importance_sampling_flag);
        }
      else
        {
          if(solPtr->importance_sampling_flag==0)
            delta=LSA_2hopt_delta_sample_estimate(solPtr,edge0fst,edge0snd,node,
                                                  edge1fst, edge1snd, realization,0);
          else
            delta=LSA_2hopt_delta_sample_estimate_window(solPtr,edge0fst,edge0snd,node,
                  edge1fst, edge1snd, realization,solPtr->importance_sampling_flag);
        }

      sum_delta +=delta;
      sum_delta_sqr+=pow(delta,2);

      if(index_realization > solPtr->minimum_realizations-1)
        {
          mean_delta=sum_delta/(double)(index_realization+1);
          sd_delta=sqrt((sum_delta_sqr
                         +(index_realization+1)*pow(mean_delta,2)
                         -2*mean_delta*sum_delta)
                        /(double)(index_realization));

          sd_estimator_delta=sd_delta/sqrt(index_realization+1);


          chebyshev_k=percentage_points_t_distribution[121][alpha_index];

          if(fabs(mean_delta) >= ((sd_estimator_delta * chebyshev_k)))
            {
              break;
            }
        }
    }


  if (mean_delta < 0.0 )
    solPtr->move_status=TRUE;

#ifdef LS_EXTRA_STATS_OUTPUT

  solPtr->solutions_explored++;
#endif

  //printf("mean_delta=%f\n",mean_delta);

  return (((double)sum_delta)/(double)GSL_MIN(index_realization+1,solPtr->maximum_realizations));
}


double
LSA_2opt_delta_sample_estimate_window(LS_List *solPtr,
                                      int edge0fst, int edge0snd,
                                      int edge1fst, int edge1snd,int index_realization, int importance_sampling)
{
  int cover0fst, cover0snd, cover1fst, cover1snd;
  long int delta = 0;
  int k,j;
  struct LS_city *this;
  const LS_DISTANCE **D = solPtr->distances;
  int window_size=solPtr->window_size;
  int no_nodes_inside_window;
  int window_length=0;

  int imp_flag;
  double correction=1.0;

  for (j=index_realization; j<index_realization+1; j++)
    {
      cover0fst = edge0fst;
      cover0snd = edge0snd;

      /* Search to the right of insertion point 0 */
      this = &solPtr->array[cover0snd];
      window_length=LSA_segment_length(solPtr,this->city,edge1snd);
      no_nodes_inside_window=(int)((window_size*solPtr->no_nodes_inside_window_percentage/100.0)+0.5);
      //printf("no_nodes_inside_window=%d\n",no_nodes_inside_window);


      while (this->city != edge1snd)
        {
          k=0;
          if((window_length < window_size) && (k < no_nodes_inside_window))
            {
              if (this->two_opt_biased_realizations[j])
                {
                  correction*=this->correction_two_opt[1];
                  break;
                }
              else
                {
                  correction*=this->correction_two_opt[0];
                  this = this->next;
                }
              k++;
            }
          else
            {
              if (this->realizations[j])
                break;
              this = this->next;
            }
#ifdef LS_EXTRA_STATS_OUTPUT
          solPtr->two_opt_scans_made++;
#endif

        }


      if (this->city == edge1snd)
        return 0.0;

      cover0snd = this->city;

      /* Search to the left of insertion point 0 */
      this = &solPtr->array[cover0fst];
      window_length=LSA_segment_length(solPtr,edge1fst,this->city);
      no_nodes_inside_window=(int)((window_size*solPtr->no_nodes_inside_window_percentage/100.0)+0.5);

      while (this->city != edge1fst)
        {
          k=0;
          if((window_length < window_size) && (k < no_nodes_inside_window))
            {
              if (this->two_opt_biased_realizations[j])
                {
                  correction*=this->correction_two_opt[1];
                  break;
                }
              else
                {
                  correction*=this->correction_two_opt[0];
                  this = this->prev;
                }
              k++;
            }
          else
            {
              if (this->realizations[j])
                break;
              this = this->prev;
            }
#ifdef LS_EXTRA_STATS_OUTPUT
          solPtr->two_opt_scans_made++;
#endif

        }

      if (this->city == edge1fst)
        return 0.0;

      cover0fst = this->city;

      cover1fst = edge1fst;
      cover1snd = edge1snd;

      /* Search to the right of insertion point 1 */
      this = &solPtr->array[cover1snd];
      window_length=LSA_segment_length(solPtr,this->city,edge0snd);
      no_nodes_inside_window=(int)((window_size*solPtr->no_nodes_inside_window_percentage/100.0)+0.5);

      while (this->city != edge0snd)
        {
          k=0;
          if((window_length < window_size) && (k < no_nodes_inside_window))
            {
              if (this->two_opt_biased_realizations[j])
                {
                  correction*=this->correction_two_opt[1];
                  break;
                }
              else
                {
                  correction*=this->correction_two_opt[0];
                  this = this->next;
                }
              k++;
            }
          else
            {
              imp_flag=0;
              if (this->realizations[j])
                break;
              this = this->next;
            }
#ifdef LS_EXTRA_STATS_OUTPUT
          solPtr->two_opt_scans_made++;
#endif

        }


      if (this->city == edge0snd)
        return 0.0;
      cover1snd = this->city;

      /* Search to the left of insertion point 1 */
      this = &solPtr->array[cover1fst];
      window_length=LSA_segment_length(solPtr,edge0fst,this->city);
      no_nodes_inside_window=(int)((window_size*solPtr->no_nodes_inside_window_percentage/100.0)+0.5);

      while (this->city != edge0fst)
        {

          k=0;
          if((window_length < window_size) && (k < no_nodes_inside_window))
            {
              if (this->two_opt_biased_realizations[j])
                {
                  correction*=this->correction_two_opt[1];
                  break;
                }
              else
                {
                  correction*=this->correction_two_opt[0];
                  this = this->prev;
                }

              k++;
            }
          else
            {
              if (this->realizations[j])
                break;
              this = this->prev;
            }

#ifdef LS_EXTRA_STATS_OUTPUT
          solPtr->two_opt_scans_made++;
#endif

        }

      if (this->city == edge0fst)
        return 0.0;
      cover1fst = this->city;

      delta =-D[cover0fst][cover0snd]
             -D[cover1fst][cover1snd]
             +D[cover0fst][cover1fst]
             +D[cover1snd][cover0snd];
    }


  return ((double)delta*correction);


}



double
LSA_2hopt_delta_sample_estimate_window(LS_List *solPtr,
                                       int edge0fst, int edge0snd, int node,
                                       int edge1fst, int edge1snd, int index_realization,int importance_sampling)
{
  int cover0fst, cover0snd, cover1fst, cover1snd;
  long int sum_delta = 0;
  int j;
  struct LS_city *this=NULL;
  const LS_DISTANCE **D = solPtr->distances;

  //	int window_size=solPtr->window_size;
  double correction=1.0;


  for (j=index_realization; j<index_realization+1; j++)
    {
      if(importance_sampling)
        {
          if (!(solPtr->array[node].two_h_opt_biased_realizations[j]))
            {

              return 0.0;
            }

          correction*=solPtr->array[node].correction_two_h_opt;
        }
      else
        {
          if (!solPtr->array[node].realizations[j])
            {
              return 0.0;
            }
        }


      cover0fst = edge0fst;
      cover0snd = edge0snd;

      /* Search to the right of insertion point 0 */
      this = &solPtr->array[cover0snd];

      while (this->city != edge1snd)
        {

          if (this->realizations[j])
            break;
          this = this->next;
#ifdef LS_EXTRA_STATS_OUTPUT

          solPtr->two_opt_scans_made++;
#endif

        }

      if (this->city == edge1snd)
        return 0.0;

      cover0snd = this->city;

      /* Search to the left of insertion point 0 */
      this = &solPtr->array[cover0fst];

      while (this->city != edge1fst)
        {
          if (this->realizations[j])
            break;
          this = this->prev;
#ifdef LS_EXTRA_STATS_OUTPUT

          solPtr->two_opt_scans_made++;
#endif

        }


      if (this->city == edge1fst)
        return 0.0;

      cover0fst = this->city;

      cover1fst = edge1fst;
      cover1snd = edge1snd;


      /* Search to the right of insertion point 1 */
      this = &solPtr->array[cover1snd];

      while (this->city != edge0snd)
        {
          if (this->realizations[j])
            break;
          this = this->next;
#ifdef LS_EXTRA_STATS_OUTPUT

          solPtr->two_opt_scans_made++;
#endif

        }

      if (this->city == edge0snd)
        return 0.0;
      cover1snd = this->city;

      /* Search to the left of insertion point 1 */
      this = &solPtr->array[cover1fst];

      while (this->city != edge0fst)
        {

          if (this->realizations[j])
            break;
          this = this->prev;
#ifdef LS_EXTRA_STATS_OUTPUT

          solPtr->two_opt_scans_made++;
#endif

        }

      if (this->city == edge0fst)
        return 0.0;
      cover1fst = this->city;


      sum_delta +=-D[cover0fst][cover0snd]
                  -D[cover1fst][node]
                  -D[node][cover1snd]
                  +D[cover0fst][node]
                  +D[node][cover0snd]
                  +D[cover1fst][cover1snd];
    }


  return ((double)sum_delta*correction);

}












