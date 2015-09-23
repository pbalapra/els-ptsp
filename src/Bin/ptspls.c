/*
#####################################################################
##########    Estimation-based Local Search for the PTSP    #########
#####################################################################

      Version: 1.0
      File:    ptspls.c
      Authors: Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
      Purpose: Main functions such as command line parser, initialization and start  
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


#include <error.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <argp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <unistd.h>

#include "heuristics.h"
#include "readFile.h"
#include "sampleLS.h"
#include "stopwatch.h"

#ifdef TESTLS_SOLUTION_LONGINT
#  define TESTLS_SOLUTION_INT long int
#else
#  define TESTLS_SOLUTION_INT int
#endif

#if defined(TESTLS_DISTANCE_LONGINT) && defined(TESTLS_DISTANCE_DOUBLE)
#error "-DTESTLS_DISTANCE_LONGINT conflicts with -DTESTLS_DISTANCE_DOUBLE"
#endif

#ifdef TESTLS_DISTANCE_LONGINT
#  define TESTLS_DISTANCE long int
#elif TESTLS_DISTANCE_DOUBLE
#  define TESTLS_DISTANCE double
#else
#  define TESTLS_DISTANCE int
#endif



typedef enum {LS_ESTIMATE,
              LS_APPROXIMATE} LS_APPROACH_TYPE;

typedef enum {LS_NN_EXPLORATION,
              LS_QNN_EXPLORATION} LS_EXPLORATION_TYPE;

typedef enum {LS_HOMOGENEOUS,
              LS_HETEROGENEOUS} LS_PROBLEM_TYPE;

typedef enum {LS_INIT_NI,
              LS_INIT_FI,
              LS_INIT_NN,
              LS_INIT_RS,
              LS_INIT_SF} LS_INIT_TYPE;

#define LS_MAX_VERBOSE 3
#define LS_MAX_VERBOSE_STRING "3"

#define LS_MAX_SEED ULONG_MAX
#define LS_MAX_REALIZATIONS 10000
#define LS_MAX_NEIGHBORS 1000
#define LS_MAX_TIME 3153600.0
#define LS_MAX_ITERATIONS INT_MAX

#define LS_DEFAULT_SEED 0
#define LS_DEFAULT_REALIZATIONS 5
#define LS_DEFAULT_DEPTH 0
#define LS_DEFAULT_SAMPLING_TYPE 1
#define LS_DEFAULT_NEIGHBORS 40
#define LS_DEFAULT_TIME 0.0000
#define LS_DEFAULT_ITERATIONS 0
#define LS_DEFAULT_ALPHA 0.05
#define LS_DEFAULT_IS 1
#define LS_DEFAULT_DELTAPROB 0.069
#define LS_DEFAULT_DELTADASHPROB 0.57
#define LS_DEFAULT_WINDOWSIZE 1.3
#define LS_DEFAULT_NODES 10

#define LSOPTION_APPROACH_ESTIMATION "estimation"
#define LSOPTION_APPROACH_APPROXIMATION "approximation"

#define LSOPTION_ESTIMATION_FIXED "fixed-estimation"
#define LSOPTION_ESTIMATION_ADAPTIVE_TYPE_1 "adaptive sampling by t-test"

#define LSOPTION_ESTIMATION_IMPSAMP_TYPE_0 "do not use importance sampling"
#define LSOPTION_ESTIMATION_IMPSAMP_TYPE_1 "biasing only the nodes close to the exchange moves with respect to window size"



#define LSOPTION_EXPLORATION_QNN "quadrant-nearest-neighbor"
#define LSOPTION_EXPLORATION_NN "nearest-neighbor"


#define LSOPTION_PROBLEM_HOMOGENEOUS "homogeneous"
#define LSOPTION_PROBLEM_HETEROGENEOUS "heterogeneous"

#define LSOPTION_INIT_NI "nearest-insertion"
#define LSOPTION_INIT_FI "farthest-insertion"
#define LSOPTION_INIT_NN "nearest-neighbor"
#define LSOPTION_INIT_RS "radial-sort"
#define LSOPTION_INIT_SF "space-filling"

#define LSOPTION_INIT_SHORT_NI "NI"
#define LSOPTION_INIT_SHORT_FI "FI"
#define LSOPTION_INIT_SHORT_NN "NN"
#define LSOPTION_INIT_SHORT_RS "RS"
#define LSOPTION_INIT_SHORT_SF "SF"


#define LS_LONGOPTIONONLY_ITERATIONS 1
#define LS_LONGOPTIONONLY_INITIALIZATION 2
#define LS_LONGOPTIONONLY_NEIGHBORHOOD 3

#define LS_LONGOPTIONONLY_ESTIMATION_IMP_SAMPLING 4
#define LS_LONGOPTIONONLY_ESTIMATION_DELTA_PROB 5
#define LS_LONGOPTIONONLY_ESTIMATION_DELTA_DASH_PROB 6
#define LS_LONGOPTIONONLY_ESTIMATION_WINDOW_SIZE 7
#define LS_LONGOPTIONONLY_ESTIMATION_WINDOW_SIZE_NODES 8


#define EITHER(a,b) "Either '" a  "' [default] or '" b "'"


typedef struct
  {
    int is_approach_estimation;
    int is_approach_approximation;
    int is_exploration_qnn;
    int is_exploration_nn;
    int is_exploration_exaustive;
    int is_approach_importance_sampling;
  }
LS_CHECK_TYPE;

LS_CHECK_TYPE check = {0,0,0,0};

const char *argp_program_version = "ptspls 0.0";
const char *argp_program_bug_address = "<prasanna@iridia.ulb.ac.be,prasannaprakash@gmail.com>";

/* Program documentation. */
static char doc[] = "ptspls -- Local search for the PTSP";
static char args_doc[] = "INSTANCE";

static struct argp_option options[] =
    {
      {"verbose",
        'v',
        "L",
        OPTION_ARG_OPTIONAL,
        "Verbosity level from 0 (silent) [default] to "
        LS_MAX_VERBOSE_STRING " (excessive).  Without the argument, "
        "the current verbose level is increased of a unit"
      },
      {"initialization",
       LS_LONGOPTIONONLY_INITIALIZATION,
       "I",
       0,
       "Initialization heuristic.  One of "
       "'" LSOPTION_INIT_SHORT_NN "' (" LSOPTION_INIT_NN ") [default], "
       "'" LSOPTION_INIT_SHORT_NI "' (" LSOPTION_INIT_NI "), "
       "'" LSOPTION_INIT_SHORT_FI "' (" LSOPTION_INIT_FI "), "
       "'" LSOPTION_INIT_SHORT_RS "' (" LSOPTION_INIT_RS "), or "
       "'" LSOPTION_INIT_SHORT_SF "' (" LSOPTION_INIT_SF ")" },
      {"exploration",
       'e',
       "E",
       0,
       "The method to be used for exploring the neighborhood. One of "
       "'" LSOPTION_EXPLORATION_QNN "' (" LSOPTION_EXPLORATION_QNN ") [default], "
       "'" LSOPTION_EXPLORATION_NN "' (" LSOPTION_EXPLORATION_NN "), "
      },
      {"problem",
       'p',
       "P",
       0,
       "The kind of PTSP problem instance to be solved.  "
       EITHER(LSOPTION_PROBLEM_HETEROGENEOUS,LSOPTION_PROBLEM_HOMOGENEOUS)},
      {"seed",
       's',
       "S",
       0,
       "Seed for initializing the random number generator"},

      {"sampling",
       'k',
       "K",
       0,
       "The type of sampling to be used in the estimation approach: One of "
       "'0' (" LSOPTION_ESTIMATION_FIXED ") , or "
       "'1' (" LSOPTION_ESTIMATION_ADAPTIVE_TYPE_1 ") [default] " },
      {"importance-sampling",
       LS_LONGOPTIONONLY_ESTIMATION_IMP_SAMPLING,
       "1",
       0,
       "\nThe type of importance sampling to be used in the estimation approach: One of "
       "'0' (" LSOPTION_ESTIMATION_IMPSAMP_TYPE_0 ") , or "
       "'1' (" LSOPTION_ESTIMATION_IMPSAMP_TYPE_1 ") [default] " },
      {"realizations",
       'r',
       "R",
       0,
       "Number of realizations: "
       "For k=0, this is the fixed number of realizations and "
       "For k>0, this is the minimum number of realizations"},
      {0,
       0,
       0,
       0,
       "If --sampling=1" },
      {"alpha",
       'c',
       "C",
       0,
       "The value of alpha used in the t-test (significance level). Choose one of : 0.01, 0.02, 0.05, 0.10 "
      },
      {0,
       0,
       0,
       0,
       "If --importance-sampling=1" },
      {"deltaProbability",
       LS_LONGOPTIONONLY_ESTIMATION_DELTA_PROB,
       "deltaP",
       0,
       "\nThe biased value of probabiltiy for 2-exchange moves."
      },
      {"deltaDashProbability",
       LS_LONGOPTIONONLY_ESTIMATION_DELTA_DASH_PROB,
       "deltaDashP",
       0,
       "The biased value of probabiltiy for node insertion moves."
      },
      {"windowsize",
       LS_LONGOPTIONONLY_ESTIMATION_WINDOW_SIZE,
       "windowSize",
       0,
       "\nThe size of the window in 2-exchange moves for importance sampling in percentage [0 to 20]."
      },
      {"nodes",
       LS_LONGOPTIONONLY_ESTIMATION_WINDOW_SIZE_NODES,
       "nodes",
       0,
       "\nThe number of nodes inside the window in percentage [0 to 50]."
      },
      {"neighbors",
       'n',
       "N",
       0,
       "Number of neighbors in the nearest neighbor or quadrant-nearest-neighbor exploration. In the case of quadrant-nearest-neighbor this number must be a multiple of 4."},
      {0,
       0,
       0,
       0,
       "Stopping criterion"},
      {"time",
       't',
       "T",
       0.000000,
       "Stop after T seconds"},
      {"steps",
       LS_LONGOPTIONONLY_ITERATIONS,
       "N",
       0,
       "Stop after N steps"},
      {0}
    };





struct arguments
  {
    int verbose;
    FILE* instance_file;
    char* instance_file_name;
    LS_INIT_TYPE initialization;
    LS_APPROACH_TYPE approach;
    LS_EXPLORATION_TYPE exploration;
    LS_PROBLEM_TYPE problem;
    unsigned long int seed;
    int realizations;
    int sampling_type;
    int neighbors;
    double time;
    int iterations;
    float alpha;
    int importance_sampling;
    float deltaProb;
    float deltaDashProb;
    float windowsize;
    float nodes;
  };

static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the INPUT argument from `argp_parse', which we
     know is a pointer to our arguments structure. */
  struct arguments *arguments = (struct arguments*)state->input;

  switch (key)
    {
    case 'v':
      if (arg)
        arguments->verbose = atoi(arg);
      else
        arguments->verbose++;
      if (arguments->verbose<0 || arguments->verbose>LS_MAX_VERBOSE)
        argp_error(state,"Illegal verbose level");
      break;

    case LS_LONGOPTIONONLY_INITIALIZATION:
      if (strcmp(arg,LSOPTION_INIT_NN)==0 ||
          strcmp(arg,LSOPTION_INIT_SHORT_NN)==0)
        arguments->initialization = LS_INIT_NN;
      else if (strcmp(arg,LSOPTION_INIT_NI)==0 ||
               strcmp(arg,LSOPTION_INIT_SHORT_NI)==0)
        arguments->initialization = LS_INIT_NI;
      else if (strcmp(arg,LSOPTION_INIT_FI)==0 ||
               strcmp(arg,LSOPTION_INIT_SHORT_FI)==0)
        arguments->initialization = LS_INIT_FI;
      else if (strcmp(arg,LSOPTION_INIT_RS)==0 ||
               strcmp(arg,LSOPTION_INIT_SHORT_RS)==0)
        arguments->initialization = LS_INIT_RS;
      else if (strcmp(arg,LSOPTION_INIT_SF)==0 ||
               strcmp(arg,LSOPTION_INIT_SHORT_SF)==0)
        arguments->initialization = LS_INIT_SF;
      else
        argp_error(state,
                   "Illegal argument '%s' for option '%s'",
                   arg, "initialization");
      break;

    case 'a':
      if (strcmp(arg,LSOPTION_APPROACH_ESTIMATION)==0)
        {
          if (check.is_approach_approximation)
            argp_error(state,"Options are not consistent");
          check.is_approach_estimation=1;
          arguments->approach = LS_ESTIMATE;
        }
      else if (strcmp(arg,LSOPTION_APPROACH_APPROXIMATION)==0)
        {
          if (check.is_approach_estimation)
            argp_error(state,"Options are not consistent");
          check.is_approach_approximation=1;
          arguments->approach = LS_APPROXIMATE;
        }
      else
        argp_error(state,
                   "Illegal argument '%s' for option '%s'",
                   arg, "approach");
      break;



    case 'e':
      if (strcmp(arg,LSOPTION_EXPLORATION_NN)==0)
        {
          if (check.is_exploration_qnn)
            argp_error(state,"Options are not consistent");
          check.is_exploration_nn=1;
          arguments->exploration = LS_NN_EXPLORATION;
        }
      else if (strcmp(arg,LSOPTION_EXPLORATION_QNN)==0)
        {
          if (check.is_exploration_nn)
            argp_error(state,"Options are not consistent");
          check.is_exploration_qnn=1;
          arguments->exploration = LS_QNN_EXPLORATION;
        }
      else
        argp_error(state,
                   "Illegal argument '%s' for option '%s'",
                   arg, "exploration");
      break;
    case 'p':
      if (strcmp(arg,LSOPTION_PROBLEM_HOMOGENEOUS)==0)
        arguments->problem = LS_HOMOGENEOUS;
      else if (strcmp(arg,LSOPTION_PROBLEM_HETEROGENEOUS)==0)
        arguments->problem = LS_HETEROGENEOUS;
      else
        argp_error(state,
                   "Illegal argument '%s' for option '%s'",
                   arg, "problem");
      break;

    case 's':
      arguments->seed = strtoul(arg,NULL,10);
      if (atoi(arg)<0 || arguments->seed>LS_MAX_SEED)
        argp_error(state,"Illegal seed");
      break;
    case 'k':
      if (check.is_approach_approximation)
        argp_error(state,"Options are not consistent");
      check.is_approach_estimation=1;
      arguments->sampling_type = atoi(arg);
      if (arguments->sampling_type < 0 ||
          arguments->sampling_type > 3)
        argp_error(state,"Illegal type of sampling");
      break;
    case 'r':
      if (check.is_approach_approximation)
        argp_error(state,"Options are not consistent");
      check.is_approach_estimation=1;
      arguments->realizations = atoi(arg);
      if (arguments->realizations<1 ||
          arguments->realizations>LS_MAX_REALIZATIONS)
        argp_error(state,"Illegal number of realizations");
      break;
    case 'c':
      if (check.is_approach_approximation)
        argp_error(state,"Options are not consistent");
      check.is_approach_estimation=1;
      arguments->alpha = atof(arg);
      if (arguments->alpha<0.0 ||
          arguments->alpha>10.0)
        argp_error(state,"Illegal alpha");
      break;

    case 'n':
      if (check.is_exploration_exaustive)
        argp_error(state,"Options are not consistent");
      if(check.is_exploration_nn)
        check.is_exploration_nn=1;
      else
        check.is_exploration_qnn=1;
      arguments->neighbors = atoi(arg);
      if (arguments->neighbors<1 ||
          arguments->neighbors>LS_MAX_NEIGHBORS)
        argp_error(state,"Illegal number of neighbors");
      break;

    case 't':
      arguments->time = atof(arg);
      if (arguments->time<0.00000 ||
          arguments->time>LS_MAX_TIME)
        argp_error(state,"Illegal time");
      break;

    case LS_LONGOPTIONONLY_ITERATIONS:
      arguments->iterations = atoi(arg);
      if (arguments->iterations<0 ||
          arguments->iterations>LS_MAX_ITERATIONS)
        argp_error(state,"Illegal number of iterations");
      break;

    case LS_LONGOPTIONONLY_ESTIMATION_IMP_SAMPLING:
      if (check.is_approach_approximation)
        argp_error(state,"Options are not consistent");
      check.is_approach_estimation=1;
      arguments->importance_sampling = atoi(arg);
      if (arguments->importance_sampling<0 ||
          arguments->importance_sampling>5)
        argp_error(state,"Illegal option for importance sampling");
      break;



    case LS_LONGOPTIONONLY_ESTIMATION_DELTA_PROB:
      if (check.is_approach_approximation)
        argp_error(state,"Options are not consistent");
      check.is_approach_estimation=1;
      //check.is_approach_importance_sampling>0;
      arguments->deltaProb = atof(arg);
      if (arguments->deltaProb<0.0 ||
          arguments->deltaProb>1.0)
        argp_error(state,"Illegal option for deltaProbability");
      break;

    case LS_LONGOPTIONONLY_ESTIMATION_DELTA_DASH_PROB:
      if (check.is_approach_approximation)
        argp_error(state,"Options are not consistent");
      check.is_approach_estimation=1;
      check.is_approach_importance_sampling=3;
      arguments->deltaDashProb = atof(arg);
      if (arguments->deltaDashProb<0.0 ||
          arguments->deltaDashProb>1.0)
        argp_error(state,"Illegal option for deltaDashProbability");
      break;

    case LS_LONGOPTIONONLY_ESTIMATION_WINDOW_SIZE:
      if (check.is_approach_approximation)
        argp_error(state,"Options are not consistent");
      check.is_approach_estimation=1;
      check.is_approach_importance_sampling=5;
      arguments->windowsize = atof(arg);
      if (arguments->windowsize<0 ||
          arguments->windowsize>50)
        argp_error(state,"Illegal option for window size");
      break;

    case LS_LONGOPTIONONLY_ESTIMATION_WINDOW_SIZE_NODES:
      if (check.is_approach_approximation)
        argp_error(state,"Options are not consistent");
      check.is_approach_estimation=1;
      check.is_approach_importance_sampling=5;
      arguments->nodes = atof(arg);
      if (arguments->nodes<0 ||
          arguments->nodes>100)
        argp_error(state,"Illegal option for nodes percentage");
      break;

    case ARGP_KEY_ARG:
      if (state->arg_num >=1)
        /* Too many arguments */
        argp_error(state,"Too many arguments");
      arguments->instance_file = fopen(arg,"r");
      if (!arguments->instance_file)
        argp_failure(state,EXIT_FAILURE,0,"Cannot open %s for reading",arg);
      arguments->instance_file_name = strdup(arg);
      if (!arguments->instance_file_name)
        argp_failure(state,EXIT_FAILURE,0,
                     "Cannot allocate memory for storing instance name");
      break;

    case ARGP_KEY_END:
      if (state->arg_num < 1)
        /* Not enough arguments */
        argp_error(state,"Instance file is missing");
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

void
check_arguments(struct arguments arguments)
{
  int no_cities=0;

  if ((check.is_exploration_nn &&
       arguments.exploration!=LS_NN_EXPLORATION) ||
      (check.is_exploration_qnn &&
       arguments.exploration!=LS_QNN_EXPLORATION) ||
      (check.is_approach_approximation &&
       arguments.approach!=LS_APPROXIMATE) ||
      (check.is_approach_estimation &&
       arguments.approach!=LS_ESTIMATE))
    error(EXIT_FAILURE,0,"Options are not consistent");

  if(arguments.time<=0.0)
    arguments.time=1.0; /*to be removed*/

  if (!(arguments.time) && !arguments.iterations)
    error(EXIT_FAILURE,0,"No stopping criterion is given");

  /* read the number of cities from instance_file */
  while (!feof(arguments.instance_file))
    if (fscanf(arguments.instance_file,"DIMENSION : %d",&no_cities))
      break;
    else
      fscanf(arguments.instance_file,"%*[^\n]\n");
  if (feof(arguments.instance_file))
    error(EXIT_FAILURE,0,
          "Error parsing instance file %s",arguments.instance_file_name);
  rewind(arguments.instance_file);



  if(check.is_exploration_qnn)
    if (arguments.neighbors%4!=0)
      error(EXIT_FAILURE,0,"Illegal number of neighbours for quadrant-nearest-neighbor. It must be a multiple of 4. ");

  if (arguments.neighbors>no_cities-1)
    error(EXIT_FAILURE,0,"Illegal number of neighbors");

  if (arguments.sampling_type>1)
    error(EXIT_FAILURE,0,"Illegal sampling type for estimation");

  if (arguments.alpha>20.0)
    error(EXIT_FAILURE,0,"Illegal alpha for adaptive estimation");

  if(arguments.importance_sampling > 1)
    error(EXIT_FAILURE,0,"Illegal type for importance sampling");

  if(arguments.importance_sampling == 1)
    if (arguments.deltaProb<0.0 || arguments.deltaProb>1.0  )
      error(EXIT_FAILURE,0,"Illegal deltaProbability for importance sampling");

  if(arguments.importance_sampling ==1)
    if (arguments.deltaDashProb<0.0 || arguments.deltaDashProb>1.0  )
      error(EXIT_FAILURE,0,"Illegal deltaDashProbability for importance sampling");

}

void
print_parameters(struct arguments arguments)
{
  printf("instance: %s\n",arguments.instance_file_name);

  printf("initialization: ");
  switch (arguments.initialization)
    {
    case LS_INIT_NI:
      printf("%s\n",LSOPTION_INIT_NI);
      break;
    case LS_INIT_FI:
      printf("%s\n",LSOPTION_INIT_FI);
      break;
    case LS_INIT_NN:
      printf("%s\n",LSOPTION_INIT_NN);
      break;
    case LS_INIT_RS:
      printf("%s\n",LSOPTION_INIT_RS);
      break;
    case LS_INIT_SF:
      printf("%s\n",LSOPTION_INIT_SF);
      break;
    default:
      printf("??\n");
      error(EXIT_FAILURE,0,"Illegal initialization heuristic");
    }

  printf("approach: ");
  switch (arguments.approach)
    {
    case LS_ESTIMATE:
      printf("%s\n",LSOPTION_APPROACH_ESTIMATION);
      break;
    case LS_APPROXIMATE:
      printf("%s\n",LSOPTION_APPROACH_APPROXIMATION);
      break;
    default:
      printf("??\n");
      error(EXIT_FAILURE,0,"Illegal approach");
    }

  printf("exploration: ");
  switch (arguments.exploration)
    {
    case LS_QNN_EXPLORATION:
      printf("%s\n",LSOPTION_EXPLORATION_QNN);
      break;
    case LS_NN_EXPLORATION:
      printf("%s\n",LSOPTION_EXPLORATION_NN);
      break;
    default:
      printf("??\n");
      error(EXIT_FAILURE,0,"Illegal exploration");
    }

  printf("problem: ");
  switch (arguments.problem)
    {
    case LS_HOMOGENEOUS:
      printf("%s\n",LSOPTION_PROBLEM_HOMOGENEOUS);
      break;
    case LS_HETEROGENEOUS:
      printf("%s\n",LSOPTION_PROBLEM_HETEROGENEOUS);
      break;
    default:
      printf("??\n");
      error(EXIT_FAILURE,0,"Illegal problem type");
    }

  printf("seed: %lu\n",arguments.seed);

  if (arguments.approach==LS_ESTIMATE)
    printf("realizations: %d\n",arguments.realizations);

  if (arguments.exploration==LS_NN_EXPLORATION||arguments.exploration==LS_QNN_EXPLORATION)
    printf("neighbors: %d\n",arguments.neighbors);

  if (arguments.time)
    printf("time: %lf\n",arguments.time);

  if (arguments.iterations)
    printf("steps: %d\n",arguments.iterations);

  if (arguments.alpha)
    printf("alpha: %f\n",arguments.alpha);


  if(arguments.importance_sampling)
    printf("Importance sampling type: %d\n",arguments.importance_sampling);

  if (arguments.deltaProb)
    printf("deltaProbabiltiy: %f\n",arguments.deltaProb);

  if (arguments.deltaDashProb)
    printf("deltaDashProbabiltiy: %f\n",arguments.deltaDashProb);

  if (arguments.windowsize)
    printf("windowsize: %f\n",arguments.windowsize);

  printf("verbose: %d\n",arguments.verbose);
}


void
LS_instance_free(problem *insPtr)
{
  int i;
  for ( i= 0  ; i < insPtr->n ; i++ )
    {
      free(insPtr->distance[i]);
    }
  free(insPtr->distance);
  free(insPtr->nodeptr);
}


int
main(int argc, char **argv)
{
  struct arguments arguments =
      {
        0, /* verbose */
        0, /* instance */
        NULL, /* instance file name */
        LS_INIT_NN, /* initialization */
        LS_ESTIMATE, /* approach */
        LS_QNN_EXPLORATION, /* exploration */
        LS_HETEROGENEOUS, /* problem */
        LS_DEFAULT_SEED, /* seed */
        LS_DEFAULT_REALIZATIONS, /* realizations */
        LS_DEFAULT_SAMPLING_TYPE,
        LS_DEFAULT_NEIGHBORS, /* neighbors */
        LS_DEFAULT_TIME, /* time */
        LS_DEFAULT_ITERATIONS, /* iterations */
        LS_DEFAULT_ALPHA,
        LS_DEFAULT_IS,
        LS_DEFAULT_DELTAPROB,
        LS_DEFAULT_DELTADASHPROB,
        LS_DEFAULT_WINDOWSIZE,
        LS_DEFAULT_NODES
      };

  argp_parse (&argp, argc, argv, 0, 0, &arguments);

  check_arguments(arguments);
  stopwatch_start(arguments.time);
  //  printf("%f\n",arguments.time);


  if (arguments.verbose)
    print_parameters(arguments);

  //  print_parameters(arguments);
  //  exit(0);

  /* do the rest... */
  int no_steps = arguments.iterations;
  int no_neighbors=arguments.neighbors;
  double time = arguments.time;
  int no_realizations=arguments.realizations;
  problem instance;
  int no_cities,i;
  int sampling_type;
  double *prob_vec;
  int homoflag;

  LS_List solution;


  instance.name=arguments.instance_file_name;
  instance.nodeptr=read_ptsp(arguments.instance_file_name, &instance);
  instance.distance=compute_distances(&instance);
  no_cities=instance.n;
  sampling_type=arguments.sampling_type;





  gsl_rng *R;
  R = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set (R, (unsigned long int)arguments.seed);

  prob_vec = (double *) malloc(no_cities * sizeof(double));
  for(int city=0; city<no_cities; city++)
    {

      if(instance.nodeptr[city].prob > 0.99)
        prob_vec[city]=0.99;
      else
        prob_vec[city]=instance.nodeptr[city].prob;
    }




  switch (arguments.initialization)
    {
    case LS_INIT_NI:
      {
        /*printf("%s\n",LSOPTION_INIT_NI);*/
        initialSolution=nearestInsertion;
      }
      break;
    case LS_INIT_FI:
      {
        /*printf("%s\n",LSOPTION_INIT_FI);*/
        initialSolution=farthestInsertion;
      }
      break;
    case LS_INIT_NN:
      {
        /*printf("%s\n",LSOPTION_INIT_NN);*/
        initialSolution=nearestNeighbor;
      }
      break;
    case LS_INIT_RS:
      {
        /*printf("%s\n",LSOPTION_INIT_RS);*/
        initialSolution=radialSort;
      }
      break;
    case LS_INIT_SF:
      {
        /*printf("%s\n",LSOPTION_INIT_SF);*/
        initialSolution=spaceFilling;
      }
      break;
    default:
      printf("??\n");
      error(EXIT_FAILURE,0,"Illegal initialization heuristic");
    }


  switch (arguments.problem)
    {
    case LS_HOMOGENEOUS:
      {
        homoflag=1;
      }
      break;
    case LS_HETEROGENEOUS:
      {
        homoflag=0;
      }
      break;
    default:
      printf("??\n");
      error(EXIT_FAILURE,0,"Illegal problem type");
    }

  long int* apriori_solution=initialSolution(&instance);
  printf("Step\t%3d\t",0);
  printf("Total_Time\t%5.8f\t", 0.00001);
  printf("Tour\t");
  for (i=0; i<no_cities; i++)
    printf("%5ld",apriori_solution[i]);
  printf("\n");

  printf("Step\t%3d\t",0);
  printf("Total_Time\t%5.8f\t", 0.00001);
  printf("Tour\t");
  for (i=0; i<no_cities; i++)
    printf("%5ld",apriori_solution[i]);
  printf("\n");
  
  printf("Step\t%3d\t",0);
  printf("Total_Time\t%5.8f\t", 0.00001);
  printf("Tour\t");
  for (i=0; i<no_cities; i++)
    printf("%5ld",apriori_solution[i]);
  printf("\n");


  switch (arguments.approach)
    {
    case LS_ESTIMATE:
      {
        solution = LS_solution_allocate(no_cities,
                                        no_realizations,
                                        prob_vec,
                                        instance.distance,
                                        arguments.alpha,
                                        arguments.importance_sampling,
                                        arguments.deltaProb,
                                        arguments.deltaDashProb,
                                        arguments.windowsize,
                                        arguments.nodes
                                       );
        /* Add realizations */
        LS_resample_realizations(solution,R,sampling_type);


        LS_solution_set(&solution,apriori_solution);



        switch (arguments.exploration)
          {
          case LS_QNN_EXPLORATION:
            LS_solution_sort_quad_neighbors(&instance,&solution,no_neighbors);
            break;
          case LS_NN_EXPLORATION:
            LS_solution_sort_neighbors(&solution,no_neighbors);
            break;
          default:
            printf("??\n");
            error(EXIT_FAILURE,0,"Illegal exploration");
          }

        LS_2hnndlbfls_times(no_steps,&solution,R,time,arguments.verbose,sampling_type);


        LS_solution_free(&solution);
      }
      break;
    default:
      printf("??\n");
      error(EXIT_FAILURE,0,"Illegal approach");
    }



  gsl_rng_free(R);
  free(prob_vec);
  free(apriori_solution);
  LS_instance_free(&instance);
  free(arguments.instance_file_name);
  fclose(arguments.instance_file);
  exit(EXIT_SUCCESS);

}


