/*
 
 ######################################################
 ########## Evaluator for the PTSP solutions ##########
 ######################################################
 
      Version: 1.0
      File:    evaluate.h
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

#define LINE_BUF_LEN 100
#define TRACE( x )
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
#define TRUE 1
#define FALSE 0

struct point
  {
    double x;
    double y;
    double p;
  };

struct problem
  {
    char          name[LINE_BUF_LEN];      	 /* instance name */
    char          edge_weight_type[LINE_BUF_LEN];  /* selfexplanatory */
    long int      n;                      /* number of cities */
    struct point  *nodeptr;               /* array of structs containing coordinates of nodes */
    long int      **distance;	        	/* distance matrix: distance[i][j] gives distance
    								   between city i und j */
  };
extern struct problem instance;


long int n;
long int *currentnodeptr;
long int seed = 12345678;

long int
(*distance)(long int, long int);  /* pointer to function returning distance */

long int
round_distance (long int i, long int j);

long int
ceil_distance (long int i, long int j);

long int
geo_distance (long int i, long int j);

long int
att_distance (long int i, long int j);

struct point
      *read_ptsp(const char *ptsp_file_name);

long int
** compute_distances(void);

long int
compute_expected_cost(long int *t);

void
printTour( long int *t );


