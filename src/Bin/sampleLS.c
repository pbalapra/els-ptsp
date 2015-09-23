/*
#####################################################################
##########    Estimation-based Local Search for the PTSP    #########
#####################################################################

      Version: 1.0
      File:    sampleLS.c
      Authors: Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
      Purpose: Functions that deal with the neighborhood exploration  
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
#include <gsl/gsl_sort_int.h>
#include <search.h>

#ifdef LS_DISTANCE_LONGINT
#  include <gsl/gsl_sort_long.h>
#elif LS_DISTANCE_DOUBLE
#  include <gsl/gsl_sort_double.h>
#else
#  include <gsl/gsl_sort_int.h>
#endif

#include "sampleLS.h"
#include "stopwatch.h"
#include "adaptiveSampling.h"
#define TRUE	1
#define FALSE	0

LS_List
LS_solution_allocate_aux(int no_cities, int no_realizations,
                         const double *prob_vec, const LS_DISTANCE **D, double alpha, int importance_sampling, float deltaProb,
                         float deltaDashProb, float window_size_percent,float nodes_percent)
{
  LS_List solution;
  int i;

  solution.window_size=(int)(((window_size_percent*no_cities)/100.0)+0.5);
  solution.no_nodes_inside_window_percentage=nodes_percent;
  solution.importance_sampling_flag=importance_sampling;


  assert(no_cities>0);
  assert(no_realizations>0);



  solution.value=0.0;
  solution.no_neighbors=0;
  solution.first=NULL;
  solution.no_cities = no_cities;
  solution.no_realizations = no_realizations;
  solution.distances = D;
  //solution.maximum_realizations=GSL_MAX((no_cities<100)?100:no_cities,no_realizations);
  solution.maximum_realizations=1000; /*this is kept constant*/
  solution.minimum_realizations=5; /*this is kept constant*/

  solution.generated_realizations=0;
  solution.move_status=FALSE;
  solution.mean_avg_delta=0.0;
  solution.sum_avg_delta=0.0;

  solution.alpha=alpha;


  /*stats*/
#ifdef LS_EXTRA_STATS_OUTPUT

  solution.solutions_explored=0;
  solution.improvements_made=0;
  solution.samples_used=0;
  solution.two_opt_scans_made=0;
  solution.two_h_opt_scans_made=0;
#endif



  solution.array = malloc(no_cities*sizeof(struct LS_city));
  if (!solution.array)
    error(EXIT_FAILURE,0,"Cannot allocate memory");


  solution.position_array= malloc(no_cities*sizeof(int));
  if (!solution.position_array)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  solution.delta = malloc((solution.maximum_realizations)*sizeof(double));
  if (!solution.delta)
    error(EXIT_FAILURE,0,"Cannot allocate memory");



  for (i=0; i<no_cities; i++)
    {

      // solution.array[i].shiftProbabilityVector=NULL;
      // solution.array[i].geometric_biased_realizations=NULL;
      // solution.array[i].correction_vector=NULL;

      solution.array[i].city = i;


      if (prob_vec)
        {

          solution.array[i].probability = prob_vec[i];

          //printf("coming here\n");

          if(solution.importance_sampling_flag==0)
            {
              solution.array[i].deltaProbability = prob_vec[i];
              solution.array[i].deltaDashProbability = prob_vec[i];
            }
          else
            {
              solution.array[i].deltaProbability = GSL_MAX(deltaProb,prob_vec[i]);
              solution.array[i].deltaDashProbability = GSL_MAX(deltaDashProb,prob_vec[i]);
              solution.array[i].correction_two_opt[0]=((1-solution.array[i].probability)/(1-solution.array[i].deltaProbability));
              solution.array[i].correction_two_opt[1]=((solution.array[i].probability)/(solution.array[i].deltaProbability));
              solution.array[i].correction_two_h_opt=((solution.array[i].probability)/(solution.array[i].deltaDashProbability));
            }

        }
      else
        error(EXIT_FAILURE,0,"probability vector is not allocated");



      solution.array[i].neighbors = 0;


      solution.array[i].realizations =
        malloc((solution.maximum_realizations)*sizeof(int));
      if (!solution.array[i].realizations)
        error(EXIT_FAILURE,0,"Cannot allocate memory");

      solution.array[i].two_opt_biased_realizations =
        malloc((solution.maximum_realizations)*sizeof(int));
      if (!solution.array[i].two_opt_biased_realizations)
        error(EXIT_FAILURE,0,"Cannot allocate memory");

      solution.array[i].two_h_opt_biased_realizations =
        malloc((solution.maximum_realizations)*sizeof(int));
      if (!solution.array[i].two_h_opt_biased_realizations)
        error(EXIT_FAILURE,0,"Cannot allocate memory");


    }

  solution.realization_order = malloc(solution.maximum_realizations*sizeof(int));
  if (!solution.realization_order)
    error(EXIT_FAILURE,0,"Cannot allocate memory");
  for (i=0; i<solution.maximum_realizations; i++)
    solution.realization_order[i]=i;

  return solution;
}



void
LS_solution_set(LS_List *solPtr, const LS_SOLUTION_INT *apriori_solution)
{
  int i;
  struct LS_city *this, *prev;

  assert(apriori_solution[0]<solPtr->no_cities);
  assert(apriori_solution[0]>=0);

  solPtr->first =
    &(solPtr->array[apriori_solution[0]]);

  solPtr->position_array[apriori_solution[0]]=0;

  solPtr->first->prev =
    &(solPtr->array[apriori_solution[solPtr->no_cities-1]]);

  for (i=1,prev=solPtr->first; i<solPtr->no_cities; i++,prev=this)
    {
      assert(apriori_solution[i]<solPtr->no_cities);
      assert(apriori_solution[i]>=0);

      this = &solPtr->array[apriori_solution[i]];

      solPtr->position_array[apriori_solution[i]]=i;

      this->prev = prev;
      prev->next = this;
    }
  (solPtr->array[apriori_solution[solPtr->no_cities-1]]).next =
    solPtr->first;
}



void
LS_solution_add_realization(LS_List solution,
                            int realization_number,
                            const double *ran_num)
{
  int i;

  assert(realization_number>=0);
  //  assert(realization_number<solution.no_realizations);
  for (i=0;i<solution.no_cities;i++)
    {
      solution.array[i].realizations[realization_number] = (ran_num[i]<=solution.array[i].probability);
      //printf("%d\n",solution.array[i].realizations[realization_number]);
      if(solution.importance_sampling_flag>0)
        {
          //printf("%f\n",solution.array[i].deltaProbability);
          solution.array[i].two_opt_biased_realizations[realization_number] = (ran_num[i]<=solution.array[i].deltaProbability);
          solution.array[i].two_h_opt_biased_realizations[realization_number] = (ran_num[i]<=solution.array[i].deltaDashProbability);
        }


    }

}




void
LS_2opt_move(LS_List *solPtr, double delta,
             int edge0fst, int edge0snd,
             int edge1fst, int edge1snd)
{
  struct LS_city *prev;
  struct LS_city *this;
  int n=solPtr->no_cities;
  int h1=0, h2=0, h3=0, h4=0, help,h_pos=0;

  //   printf("edge0fst: %d,pos: %d\n",edge0fst,solPtr->position_array[edge0fst]);
  //   printf("edge0snd: %d,pos: %d\n",edge0snd,solPtr->position_array[edge0snd]);
  //   printf("edge1fst: %d,pos: %d\n",edge1fst,solPtr->position_array[edge1fst]);
  //   printf("edge1snd: %d,pos: %d\n",edge1snd,solPtr->position_array[edge1snd]);
  assert(!((solPtr->position_array[edge0fst]+1 )%n !=solPtr->position_array[edge0snd]));
  //     if((solPtr->position_array[edge0fst]+1 )%n !=solPtr->position_array[edge0snd])
  //    {
  //  	  printf("check edge0fst: %d\n",(solPtr->position_array[edge0fst]+1)%n);
  //  	  printf("check edge0snd: %d\n",solPtr->position_array[edge0snd]);
  //  	  exit(0);
  //    }
  assert(!((solPtr->position_array[edge1fst]+1)%n != solPtr->position_array[edge1snd]));
  //   if((solPtr->position_array[edge1fst]+1)%n != solPtr->position_array[edge1snd])
  //    if((solPtr->position_array[edge1fst]+1)%n != solPtr->position_array[edge1snd])
  //    {
  //  	  printf("check edge1fst: %d\n",(solPtr->position_array[edge1fst]+1)%n);
  //  	  printf("check edge1snd: %d\n",solPtr->position_array[edge1snd]);
  //  	  exit(0);
  //    }

  h1 = edge0fst;
  h2 = edge0snd;
  h3 = edge1fst;
  h4 = edge1snd;

  if ( solPtr->position_array[h3] < solPtr->position_array[h1] )
    {
      help = h1;
      h1 = h3;
      h3 = help;
      help = h2;
      h2 = h4;
      h4 = help;
    }

  if ( solPtr->position_array[h3] - solPtr->position_array[h2] < n / 2 + 1)
    {
      /* reverse inner part from pos[h2] to pos[h3] */
      this = &solPtr->array[h3];
      prev = &solPtr->array[h1];
      solPtr->array[h1].next = &solPtr->array[h3];
      h_pos=solPtr->position_array[prev->city];
      while (this->city != h2)
        {
          this->next = this->prev;
          this->prev = prev;
          solPtr->position_array[this->city]=(++h_pos)%n;
          prev = this;
          //printf("moving inner %d\n",this->city);
          //printf(" with rank -> %d\n",solPtr->position_array[this->city]);
          this = this->next;
        }
      solPtr->array[h2].prev = prev;
      solPtr->array[h2].next = &solPtr->array[h4];
      solPtr->array[h4].prev = &solPtr->array[h2];
      solPtr->position_array[h2]=(++h_pos)%n;
      //printf("moving inner %d\n",h2);
      //printf(" with rank -> %d\n",solPtr->position_array[h2]);
    }
  else
    {
      /* reverse outer part from pos[h4] to pos[h1] */
      this = &solPtr->array[h1];
      prev = &solPtr->array[h3];
      solPtr->array[h3].next = &solPtr->array[h1];
      h_pos=solPtr->position_array[prev->city];
      while (this->city != h4)
        {
          this->next = this->prev;
          this->prev = prev;
          solPtr->position_array[this->city]=(++h_pos)%n;
          prev = this;
          //printf("moving outer %d\n",this->city);
          //printf(" with rank -> %d\n",solPtr->position_array[this->city]);
          this = this->next;
        }
      solPtr->array[h4].prev = prev;
      solPtr->array[h4].next = &solPtr->array[h2];
      solPtr->array[h2].prev = &solPtr->array[h4];
      solPtr->position_array[h4]=(++h_pos)%n;
      //printf("moving outer %d\n",h4);
      //printf(" with rank -> %d\n",solPtr->position_array[h4]);

    }



  solPtr->value += delta;

}


double
LS_2nndlbfls_step(LS_List *solPtr, const int *order, gsl_rng *r, int sampling_type)
{
  int edge0fst, edge0snd, edge1fst=-1, edge1snd=-1;
  int this_city;
  const LS_DISTANCE **D = solPtr->distances;
  double delta;
  int i, j;
  LS_DISTANCE radius;
  int nn = solPtr->no_neighbors;
  int n = solPtr->no_cities;

  if(sampling_type==0)
    delta_evaluation=LS_delta;
  else if (sampling_type==1)
    delta_evaluation=LSA_delta_adaptive_sample;
  else if (sampling_type==2)
    delta_evaluation=LSA_delta_adaptive_sample;

  for (i=0; i<n; i++)
    {
      this_city = order[i];

      if (solPtr->array[this_city].dlb)
        continue;

      edge0fst = this_city;
      edge0snd = solPtr->array[this_city].next->city;
      radius = D[edge0fst][edge0snd];
      /* The loop starts from 1 because
      position 0 is the current city itself */
      for (j=1; j<nn; j++)
        {
          edge1fst = solPtr->array[edge0fst].neighbors[j];
          if (radius>D[edge0fst][edge1fst])
            {
              edge1snd = solPtr->array[edge1fst].next->city;

              if ( edge1fst==edge0snd || edge0fst==edge1snd )
                continue;

              delta=delta_evaluation(solPtr,edge0fst,edge0snd,-1,edge1fst,edge1snd,r,0);
              if (solPtr->move_status)
                {
                  LS_2opt_move(solPtr,delta,
                               edge0fst,edge0snd,edge1fst,edge1snd);
                  solPtr->array[edge0fst].dlb = 0;
                  solPtr->array[edge0snd].dlb = 0;
                  solPtr->array[edge1fst].dlb = 0;
                  solPtr->array[edge1snd].dlb = 0;
                  return(delta);
                }
            }
          else
            break;
        }


      edge0fst = solPtr->array[this_city].prev->city;
      edge0snd = this_city;
      radius = D[edge0fst][edge0snd];

      /* The loop starts from 1 because
      position 0 is the current city itself */
      for (j=1; j<nn; j++)
        {
          edge1snd = solPtr->array[edge0snd].neighbors[j];
          if (radius>D[edge0snd][edge1snd])
            {
              edge1fst=solPtr->array[edge1snd].prev->city;

              /* Seen in Thomas code but not fully understood...*/
              if ( edge1fst==edge0snd || edge0fst==edge1snd )
                continue;

              delta=delta_evaluation(solPtr,edge0fst,edge0snd,-1,edge1fst,edge1snd,r,0);

              if (solPtr->move_status)
                {
                  LS_2opt_move(solPtr,delta,
                               edge0fst,edge0snd,edge1fst,edge1snd);
                  solPtr->array[edge0fst].dlb = 0;
                  solPtr->array[edge0snd].dlb = 0;
                  solPtr->array[edge1fst].dlb = 0;
                  solPtr->array[edge1snd].dlb = 0;
                  return(delta);
                }
            }
          else
            break;
        }

      /* If we are here, its because no improving neighboring
      solution was found starting from `this_city' */
      solPtr->array[this_city].dlb = 1;
    }


  /* If we are here, its because no improving neighboring
     solution was found starting from any city. We than report that
     no improvement was made */
  return 0.0;
}



inline void
LS_2hopt_move(LS_List *solPtr, double delta,
              int edge0fst, int edge0snd, int node,
              int edge1fst, int edge1snd)
{

  struct LS_city *this;
  int h_pos=0;
  int n=solPtr->no_cities;

  //   printf("edge0fst: %d,pos: %d\n",edge0fst,solPtr->position_array[edge0fst]);
  //   printf("edge0snd: %d,pos: %d\n",edge0snd,solPtr->position_array[edge0snd]);
  //   printf("node: %d,pos: %d\n",node,solPtr->position_array[node]);
  //   printf("edge1fst: %d,pos: %d\n",edge1fst,solPtr->position_array[edge1fst]);
  //   printf("edge1snd: %d,pos: %d\n",edge1snd,solPtr->position_array[edge1snd]);

  assert(!((solPtr->position_array[edge0fst]+1 )%n !=solPtr->position_array[edge0snd]));
  //   if((solPtr->position_array[edge0fst]+1 )%n !=solPtr->position_array[edge0snd])
  //   {
  // 	  printf("check edge0fst: %d\n",(solPtr->position_array[edge0fst]+1)%n);
  // 	  printf("check edge0snd: %d\n",solPtr->position_array[edge0snd]);
  // 	  exit(0);
  //   }
  assert(!((solPtr->position_array[edge1fst]+1)%n != solPtr->position_array[node]));
  //   if((solPtr->position_array[edge1fst]+1)%n != solPtr->position_array[node])
  //   {
  // 	  printf("check edge1fst: %d\n",(solPtr->position_array[edge1fst]+1)%n);
  // 	  printf("check edge1snd: %d\n",solPtr->position_array[node]);
  // 	  exit(0);
  //   }
  assert(!((solPtr->position_array[node]+1)%n != solPtr->position_array[edge1snd]));
  //   if((solPtr->position_array[node]+1)%n != solPtr->position_array[edge1snd])
  //   {
  // 	  printf("check edge1fst: %d\n",(solPtr->position_array[edge1fst]+1)%n);
  // 	  printf("check edge1snd: %d\n",solPtr->position_array[node]);
  // 	  exit(0);
  //   }

  solPtr->array[edge0fst].next = &solPtr->array[node];
  solPtr->array[node].prev = &solPtr->array[edge0fst];
  solPtr->array[node].next = &solPtr->array[edge0snd];
  solPtr->array[edge0snd].prev = &solPtr->array[node];
  solPtr->array[edge1fst].next = &solPtr->array[edge1snd];
  solPtr->array[edge1snd].prev = &solPtr->array[edge1fst];
  solPtr->value += delta;


  if ( solPtr->position_array[edge0fst] < solPtr->position_array[edge1fst] )
    {
      h_pos=solPtr->position_array[edge0fst];
      this = &solPtr->array[node];
      /* Now update the position array */
      while ( this->city != edge1snd )
        {
          solPtr->position_array[this->city]=(++h_pos)%n;
          //printf("moving 2hopt %d\n",this->city);
          //printf(" with rank -> %d\n",solPtr->position_array[this->city]);
          this=this->next;
        }
    }
  else
    {
      h_pos=solPtr->position_array[edge0snd];
      this = &solPtr->array[node];
      /* Now update the position array */
      while ( this->city != edge1fst )
        {
          if (h_pos==0)
            h_pos=n;
          solPtr->position_array[this->city]=(--h_pos)%n;
          //assert(!(h_pos<0));
          //printf("moving 2hopt %d\n",this->city);
          //		  printf(" with rank -> %d\n",solPtr->position_array[this->city]);
          assert(!(h_pos<0));
          this=this->prev;
        }
    }

  //assert(LS_move_check_position_index(solPtr));

}


double
LS_2hnndlbfls_step(LS_List *solPtr, const int *order, gsl_rng *r, int sampling_type)
{
  int edge0fst, edge0snd, edge1fst=-1, edge1snd=-1, node=-1;
  int this_city;

  const LS_DISTANCE **D = solPtr->distances;

  double delta;
  int i, j;
  LS_DISTANCE radius;
  int nn = solPtr->no_neighbors;
  int n = solPtr->no_cities;

  if(sampling_type==0)
    delta_evaluation=LS_delta;
  else if (sampling_type==1)
    delta_evaluation=LSA_delta_adaptive_sample;
  else if (sampling_type==2)
    delta_evaluation=LSA_delta_adaptive_sample;


  for (i=0; i<n; i++)
    {
      this_city = order[i];

      if (solPtr->array[this_city].dlb)
        continue;

      edge0fst = this_city;
      edge0snd = solPtr->array[this_city].next->city;
      radius = D[edge0fst][edge0snd];

      /* The loop starts from 1 because
      position 0 is the current city itself */
      for (j=1; j<nn; j++)
        {
          edge1fst = solPtr->array[edge0fst].neighbors[j];
          if (radius>D[edge0fst][edge1fst])
            {
              edge1snd = solPtr->array[edge1fst].next->city;
              /* Seen in Thomas code but not fully understood...*/
              if ( edge1fst==edge0snd || edge0fst==edge1snd )
                continue;
              //delta=LS_2opt_delta(solPtr,edge0fst,edge0snd,edge1fst,edge1snd,r);
              delta=delta_evaluation(solPtr,edge0fst,edge0snd,-1,edge1fst,edge1snd,r,0);
              if (solPtr->move_status)
                {
                  LS_2opt_move(solPtr,delta,
                               edge0fst,edge0snd,edge1fst,edge1snd);
                  solPtr->array[edge0fst].dlb = 0;
                  solPtr->array[edge0snd].dlb = 0;
                  solPtr->array[edge1fst].dlb = 0;
                  solPtr->array[edge1snd].dlb = 0;
                  return(delta);
                }

              node = edge1fst;
              edge1fst = solPtr->array[node].prev->city;
              //            delta=LS_2hopt_delta(solPtr,edge0fst,edge0snd,node,edge1fst,edge1snd,r);
              delta=delta_evaluation(solPtr,edge0fst,edge0snd,node,edge1fst,edge1snd,r,1);
              if (solPtr->move_status)
                {
                  LS_2hopt_move(solPtr,delta,edge0fst,edge0snd,
                                node,edge1fst,edge1snd);
                  solPtr->array[node].dlb = 0;
                  solPtr->array[edge0fst].dlb = 0;
                  solPtr->array[edge0snd].dlb = 0;
                  solPtr->array[edge1fst].dlb = 0;
                  solPtr->array[edge1snd].dlb = 0;
                  return(delta);
                }
            }
          else
            break;
        }


      edge0fst = solPtr->array[this_city].prev->city;
      edge0snd = this_city;
      radius = D[edge0fst][edge0snd];

      /* The loop starts from 1 because
      position 0 is the current city itself */
      for (j=1; j<nn; j++)
        {
          edge1snd = solPtr->array[edge0snd].neighbors[j];
          if (radius>D[edge0snd][edge1snd])
            {
              edge1fst=solPtr->array[edge1snd].prev->city;

              /* Seen in Thomas code but not fully understood...*/
              if ( edge1fst==edge0snd || edge0fst==edge1snd )
                continue;

              //delta=LS_2opt_delta(solPtr,edge0fst,edge0snd,edge1fst,edge1snd,r);
              delta=delta_evaluation(solPtr,edge0fst,edge0snd,-1,edge1fst,edge1snd,r,0);
              if (solPtr->move_status)
                {
                  LS_2opt_move(solPtr,delta,
                               edge0fst,edge0snd,edge1fst,edge1snd);
                  solPtr->array[edge0fst].dlb = 0;
                  solPtr->array[edge0snd].dlb = 0;
                  solPtr->array[edge1fst].dlb = 0;
                  solPtr->array[edge1snd].dlb = 0;
                  return(delta);
                }

              node = edge1snd;
              edge1snd = solPtr->array[node].next->city;
              //            delta=LS_2hopt_delta(solPtr,edge0fst,edge0snd,node,edge1fst,edge1snd,r);
              delta=delta_evaluation(solPtr,edge0fst,edge0snd,node,edge1fst,edge1snd,r,1);
              if (solPtr->move_status)
                {
                  LS_2hopt_move(solPtr,delta,edge0fst,edge0snd,
                                node,edge1fst,edge1snd);
                  solPtr->array[node].dlb = 0;
                  solPtr->array[edge0fst].dlb = 0;
                  solPtr->array[edge0snd].dlb = 0;
                  solPtr->array[edge1fst].dlb = 0;
                  solPtr->array[edge1snd].dlb = 0;
                  return(delta);
                }
            }
          else
            break;
        }

      /* If we are here, its because no improving neighboring
      solution was found starting from `this_city' */
      solPtr->array[this_city].dlb = 1;
    }


  /* If we are here, its because no improving neighboring
     solution was found starting from any city. We than report that
     no improvement was made */
  return 0.0;
}


void
LS_reset_dlb(LS_List *solPtr)
{
  int i;
  for (i=0; i<solPtr->no_cities; i++)
    solPtr->array[i].dlb=0;
}






void
LS_solution_free(LS_List *solPtr)
{
  int i;
  for (i=0; i<solPtr->no_cities; i++)
    {
      free(solPtr->array[i].realizations);
      free(solPtr->array[i].two_opt_biased_realizations);
      free(solPtr->array[i].two_h_opt_biased_realizations);
      solPtr->array[i].realizations = NULL;
      if (solPtr->array[i].neighbors)
        {
          free(solPtr->array[i].neighbors);
          solPtr->array[i].neighbors = NULL;
        }

      //free(solPtr->array[i].geometric_biased_realizations);
      //solPtr->array[i].geometric_biased_realizations=NULL;

    }
  free(solPtr->array);
  solPtr->array = NULL;
  free(solPtr->position_array);
  free(solPtr->delta);
  free(solPtr->realization_order);
}


void
LS_solution_print_all(LS_List solution)
{
  int i, j;
  struct LS_city *this;

  printf("---------------------------------------\n");
  for (i=0,this=solution.first;i<solution.no_cities;i++,this=this->next)
    {
      printf("City: %5d\tnext: %5d\tprev: %5d\n",
             this->city,(this->next)->city,(this->prev)->city);
      printf("\tRealizations:\n\t");
      for (j=0; j<solution.maximum_realizations; j++)
        printf("%2d",this->realizations[j]);
      printf("\n");
      if (this->neighbors)
        {
          printf("\tNeighbors:\n\t");
          for (j=0; j<solution.no_neighbors; j++)
            printf("%2d",this->neighbors[j]);
          printf("\n");
        }
    }
}





void
LS_solution_print_aux(LS_List solution)
{
  int i;
  struct LS_city *this;
  this = solution.first;
  for (i=0; i<solution.no_cities; i++,this=this->next)
    printf("%5d",this->city);
}

void
LS_solution_print_aux1(LS_List solution)
{


  int no_cities=solution.no_cities;
  double expc = 0.0,partial=0.0;
  int i,r,s;
  double* prob_vec;
  double** Q;
  long int* t;
  struct LS_city *this;


  t = (long int *) malloc(no_cities*sizeof(long int));
  if (!t)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  this = solution.first;
  for (i=0; i<solution.no_cities; i++,this=this->next)
    {
      t[i]=this->city;

    }

  prob_vec = (double *) malloc(no_cities*sizeof(double));
  if (!prob_vec)
    error(EXIT_FAILURE,0,"Cannot allocate memory");
  for(i=0; i<no_cities; i++)
    {
      prob_vec[i]= solution.array[i].probability;
    }

  Q = malloc(no_cities*sizeof(double *));
  if(!Q)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  for(i=0; i < no_cities; i++)
    {
      Q[i] = malloc(no_cities*sizeof(double));
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

  for (i=0; i<no_cities; i++)
    {
      for(r=0; r< no_cities ; r++)
        {
          if(i != r )
            {
              partial = solution.distances[t[i]][t[r]]
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
  printf("cost=%lf\n", expc);
  free(t);
}



void
LS_solution_print(LS_List solution)
{
  printf("=>");
  LS_solution_print_aux(solution);
  printf(": %.3f\n",solution.value);
}

void
LS_solution_log(LS_List solution, int i)
{
  printf("Step\t%3d\t",i);
  printf("Total_Time\t%5.8f\t",stopwatch_read());
  printf("Tour\t");
  LS_solution_print_aux(solution);
  printf("\n");
}


void
LS_resample_realizations(LS_List solution, gsl_rng *r, int sampling_type)
{
  int i,j,no_realizations;
  double *ran_num;

  if(sampling_type==0)
    no_realizations=solution.no_realizations;
  else
    no_realizations=solution.maximum_realizations;


  ran_num = malloc(solution.no_cities*sizeof(double));
  if (!ran_num)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  for (i=0; i<no_realizations; i++)
    {
      for (j=0; j<solution.no_cities; j++)
        {
          ran_num[j]=gsl_rng_uniform(r);
        }
      LS_solution_add_realization(solution,i,ran_num);
    }
  solution.generated_realizations+=no_realizations;
  free(ran_num);
}



/* Static version of the function LS_solution_allocate that has been
declared as `extern inline' in the file sampleLS.h */
LS_List
LS_solution_allocate(int no_cities,int no_realizations,
                     const double *prob_vec, LS_DISTANCE **D,
                     double alpha, int importance_sampling,
                     float deltaProb, float deltaDashProb,
                     float window_size_percent,float nodes_percent)
{
  return LS_solution_allocate_aux(no_cities,no_realizations, prob_vec,
                                  (const LS_DISTANCE **) D, alpha,
                                  importance_sampling, (double) deltaProb, (double) deltaDashProb,
                                  window_size_percent,nodes_percent);
}


void
LS_solution_set_with_given_value(LS_List *solPtr, double value,
                                 const LS_SOLUTION_INT *apriori_solution)
{
  LS_solution_set(solPtr,apriori_solution);
  solPtr->value = value;
}

void
LS_Xnndlbfls_times(double (*LSstep)(LS_List*, const int*, gsl_rng *r, int sampling_type),
                   void (*LSresample)(LS_List,gsl_rng*, int),
                   int K, LS_List *solPtr, gsl_rng *r,
                   double time, int verbose, int sampling_type)
{
  int i;
  double m = -1.0;
  int *order;

  order = malloc(solPtr->no_cities*sizeof(int));
  for (i=0; i<solPtr->no_cities; i++)
    order[i]=i;
  gsl_ran_shuffle(r,order,solPtr->no_cities,sizeof(int));

  LS_reset_dlb(solPtr);

  for (i=1;  ((solPtr->move_status)||m<0) && (!K||i<K) ; i++)
    {
      if (LSresample)
        (*LSresample)(*solPtr,r,sampling_type);
      m = (*LSstep)(solPtr,order,r,sampling_type);
      solPtr->value=m;
      if ((stopwatch_read()>=0.0) && (solPtr->move_status))
        {
          //LS_solution_log(*solPtr,i);
          //LS_solution_print_aux1(*solPtr);
        }

      gsl_ran_shuffle(r,solPtr->realization_order,solPtr->maximum_realizations,sizeof(int));

    }
  LS_solution_log(*solPtr,i);
  LS_solution_log(*solPtr,i);
  //LS_solution_print_aux1(*solPtr);

  free(order);
}


void
LS_2nndlbfls(LS_List *solPtr, gsl_rng *r, double time, int verbose, int sampling_type)
{
  LS_Xnndlbfls_times(&LS_2nndlbfls_step,NULL,0,solPtr, r, time, verbose, sampling_type);
}

void
LS_2hnndlbfls(LS_List *solPtr, gsl_rng *r, double time, int verbose, int sampling_type)
{
  LS_Xnndlbfls_times(&LS_2hnndlbfls_step,NULL,0,solPtr, r, time, verbose,sampling_type);
}

void
LS_2nndlbfls_times(int K, LS_List *solPtr, gsl_rng *r,
                   double time, int verbose, int sampling_type)
{
  LS_Xnndlbfls_times(&LS_2nndlbfls_step,NULL,K,solPtr, r, time, verbose,sampling_type);
}

void
LS_2hnndlbfls_times(int K, LS_List *solPtr, gsl_rng *r,
                    double time, int verbose, int sampling_type)
{
  LS_Xnndlbfls_times(&LS_2hnndlbfls_step,NULL,K,solPtr, r, time, verbose,sampling_type);
}

void
LS_2nndlbfls_resample(LS_List *solPtr, gsl_rng *r, double time, int verbose, int sampling_type)
{
  LS_Xnndlbfls_times(&LS_2nndlbfls_step,&LS_resample_realizations,
                     0,solPtr, r, time, verbose,sampling_type);
}

void
LS_2hnndlbfls_resample(LS_List *solPtr, gsl_rng *r, double time, int verbose, int sampling_type)
{
  LS_Xnndlbfls_times(&LS_2hnndlbfls_step,&LS_resample_realizations,
                     0,solPtr, r, time, verbose,sampling_type);
}

void
LS_2nndlbfls_times_resample(int K, LS_List *solPtr, gsl_rng *r,
                            double time, int verbose, int sampling_type)
{
  LS_Xnndlbfls_times(&LS_2nndlbfls_step,&LS_resample_realizations,
                     K,solPtr, r, time, verbose,sampling_type);
}

void
LS_2hnndlbfls_times_resample(int K, LS_List *solPtr, gsl_rng *r,
                             double time, int verbose, int sampling_type)
{
  LS_Xnndlbfls_times(&LS_2hnndlbfls_step,&LS_resample_realizations,
                     K,solPtr, r, time, verbose,sampling_type);
}


double
LS_solution_compute_and_set_value(LS_List *solPtr)
{
  solPtr->value = LS_solution_return_value(*solPtr);
  return(solPtr->value);
}





void
LS_solution_sort_quad_neighbors(problem *insPtr, LS_List *solPtr,int nn)
{
  int i;
  solPtr->no_neighbors = nn;

  for (i=0; i<solPtr->no_cities; i++)
    solPtr->array[i].neighbors =
      LS_allocate_sort_quad_neighbors(insPtr,solPtr->distances,i,solPtr->no_cities,nn);
}

int*
LS_allocate_sort_quad_neighbors(problem *insPtr,const LS_DISTANCE **d,int node,int no_cities,int nn)
{

  int i;
  int no_in_quadrant[5];
  double n_xx, n_yy, c_xx, c_yy;
  int check_nn;
  int *distance_vector;
  int *help_vector;
  int *chosen;
  int ncities=no_cities;
  int nn_ls=nn;
  int nn_quad;
  int *m_nnear;
  /*    printf("\n computing nearest neighbor lists, "); */


  nn_ls = GSL_MIN(nn_ls,ncities);
  nn_quad=nn_ls/4;

  m_nnear = malloc(nn_ls*sizeof(int));
  if (!m_nnear)
    error(EXIT_FAILURE,0,"Cannot allocate memory");

  distance_vector = malloc(ncities*sizeof(int));
  help_vector = malloc(ncities*sizeof(int));
  chosen = malloc(ncities*sizeof(int));

  /*    for ( node = 0 ; node < ncities ; node++ ) { */ /* compute cnd-sets for all node */

  /*	m_nnear[node] = (LS_SOLUTION_INT *)(m_nnear + ncities) + node * nn_ls;*/


  for ( i = 0 ; i < ncities ; i++ )
    {  /* Copy distances from nodes to the others */
      distance_vector[i] = d[node][i];
      help_vector[i] = i;
      chosen[i] = FALSE;
    }

  distance_vector[node] = LONG_MAX;  /* city is not nearest neighbour */
  LS_sort(distance_vector, help_vector, 0, ncities-1);
  if(node<ncities-2)
    {
      n_xx = insPtr->nodeptr[node+1].x;
      n_yy = insPtr->nodeptr[node+1].y;
    }
  else
    {
      n_xx = insPtr->nodeptr[0].x;
      n_yy = insPtr->nodeptr[0].y;
    }
  no_in_quadrant[1] = 0;
  no_in_quadrant[2] = 0;
  no_in_quadrant[3] = 0;
  no_in_quadrant[4] = 0;
  check_nn = 0;
  for ( i = 0 ; i < ncities ; i++ )
    {
      if((help_vector[i]+1)==ncities)
        {
          c_xx = insPtr->nodeptr[0].x;
          c_yy = insPtr->nodeptr[0].y;
        }
      else
        {
          c_xx = insPtr->nodeptr[help_vector[i]+1].x;
          c_yy = insPtr->nodeptr[help_vector[i]+1].y;
        }
      if ( c_xx > n_xx && c_yy >= n_yy )
        {
          /* 1. Quadrant */
          if ( no_in_quadrant[1] < nn_quad )
            {
              no_in_quadrant[1] += 1;
              chosen[i] = TRUE;
              check_nn++;
            }
        }
      else if ( c_xx >= n_xx && c_yy < n_yy )
        {
          /* 2. Quadrant */
          if ( no_in_quadrant[2] < nn_quad )
            {
              no_in_quadrant[2] += 1;
              chosen[i] = TRUE;
              check_nn++;
            }
        }
      else if ( c_xx < n_xx && c_yy <= n_yy )
        {
          /* 3. Quadrant */
          if ( no_in_quadrant[3] < nn_quad )
            {
              no_in_quadrant[3] += 1;
              chosen[i] = TRUE;
              check_nn++;
            }
        }
      else if ( c_xx <= n_xx && c_yy > n_yy )
        {
          /* 4. Quadrant */
          if ( no_in_quadrant[4] < nn_quad )
            {
              no_in_quadrant[4] += 1;
              chosen[i] = TRUE;
              check_nn++;
            }
        }
    }
  if ( check_nn < 4 * nn_quad )
    {
      for ( i = 0 ; i < ncities ; i++ )
        {
          if( chosen[i] == FALSE )
            {
              chosen[i] = TRUE;
              /*  		    printf(" %d ",i); */
              check_nn++;
            }
          if ( check_nn == 4 * nn_quad )
            break;
        }
    }

  if ( check_nn != 4 * nn_quad )
    exit(0);

  check_nn = 0;
  i = 0;
  while ( i < ncities )
    {
      if ( chosen[i] )
        {
          m_nnear[check_nn] = help_vector[i];
          check_nn += 1;
          /*    		printf(" %d ",i); */
        }

      i++;
      if ( check_nn == 4 * nn_quad )
        break;
    }
  /*   	printf("\n"); */
  /*  }*/
  free(distance_vector);
  free(help_vector);
  free(chosen);
  /*    printf("\n    .. done\n"); */
  return m_nnear;

}




void LS_swap(int v[], int v2[], int i, int j)
{
  long int tmp;

  tmp = v[i];
  v[i] = v[j];
  v[j] = tmp;
  tmp = v2[i];
  v2[i] = v2[j];
  v2[j] = tmp;

}

void LS_sort(int v[], int v2[], int left, int right)
{
  long int k, last;

  if (left >= right)
    return;
  LS_swap(v, v2, left, (left + right)/2);
  last = left;
  for (k=left+1; k <= right; k++)
    if (v[k] < v[left])
      LS_swap(v, v2, ++last, k);
  LS_swap(v, v2, left, last);
  LS_sort(v, v2, left, last);
  LS_sort(v, v2, last+1, right);
}







int
LS_move_check_position_index(LS_List *solPtr)
{
  int no_error_flag=TRUE;
  int n=solPtr->no_cities;
  int i=0;
  int sum_pos=0;
  int sum_node=0;
  int first_node_pos;
  struct LS_city *this;
  int *city_index = malloc(n*sizeof(int));
  if (!city_index)
    error(EXIT_FAILURE,0,"Cannot allocate memory");
  gsl_sort_int_smallest_index((size_t*)city_index,n,solPtr->position_array,1,n);

  first_node_pos=city_index[i];
  this=&solPtr->array[first_node_pos];
  for ( i = 0 ; i < n ; i++ )
    {
      sum_pos+=city_index[i];
      sum_node+=this->city;
      //printf("%d = %d\n",city_index[i],this->city);
      if(this->city != city_index[i])
        {
          no_error_flag=FALSE;
          break;
        }
      this=this->next;
    }

  if ( sum_node != (n-1) * n / 2 )
    {
      printf(" tour must be flawed !!\n");
      no_error_flag=FALSE;
    }

  if ( sum_pos != (n-1) * n / 2 )
    {
      printf("position array must be flawed !!\n");
      no_error_flag=FALSE;
    }
  free(city_index);
  return no_error_flag;
}


void
LS_solution_sort_neighbors(LS_List *solPtr, int nn)
{
  int i;
  const LS_DISTANCE **D = solPtr->distances;
  int n = solPtr->no_cities;
  nn = GSL_MIN(n,nn);
  solPtr->no_neighbors = nn;

  for (i=0; i<solPtr->no_cities; i++)
    solPtr->array[i].neighbors = LS_allocate_sort_neighbors(D[i],n,nn);
}
int*
LS_allocate_sort_neighbors(const LS_DISTANCE *d, int no_cities, int nn)
{
  int *neighbours = malloc(nn*sizeof(LS_DISTANCE));
  if (!neighbours)
    error(EXIT_FAILURE,0,"Cannot allocate memory");
  assert(nn<=no_cities);
#ifdef LS_DISTANCE_LONGINT

  gsl_sort_long_smallest_index((size_t*)neighbours,nn,d,1,no_cities);
#elif LS_DISTANCE_DOUBLE

  gsl_sort_smallest_index((size_t*)neighbours,nn,d,1,no_cities);
#else

  gsl_sort_int_smallest_index((size_t*)neighbours,nn,d,1,no_cities);
#endif

  return neighbours;
}


double
LS_solution_get(LS_List solution, LS_SOLUTION_INT *apriori_solution)
{
  int i;
  struct LS_city *this;

  for (i=0,this=solution.first; i<solution.no_cities; i++,this=this->next)
    apriori_solution[i] = this->city;

  /* See sampleLS.h for an explanation of next line: */
#ifdef LS_REPEAT_FIRST_CITY_AT_END
  
  apriori_solution[solution.no_cities] = apriori_solution[0];
#endif

  return (solution.value);
}

double
LS_solution_return_value(LS_List solution)
                  {
                    struct LS_city *this;
                    int i, j;
                    int start, pos;
                    long int tmp;
                    double value = 0.0;
                    const LS_DISTANCE **D = solution.distances;

                    for (j=0; j<solution.no_realizations; j++)
                      {
                        tmp = 0;
                        pos = start = -1;
                        this = solution.first;
                        for (i=0; i<solution.no_cities; i++)
                          {
                            if (this->realizations[j])
                              {
                                if (start == -1)
                                  start = this->city;
                                else
                                  tmp += D[pos][this->city];
                                pos = this->city;
                              }
                            this = this->next;
                          }
                        if (start != -1)
                          tmp += D[pos][start];
                        value += ((double)tmp)/((double)solution.no_realizations);
                      }
                    return(value);
                  }
