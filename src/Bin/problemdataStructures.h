/*
#####################################################################
##########    Estimation-based Local Search for the PTSP    #########
#####################################################################

      Version: 1.0
      File:    problemdataStructures.h
      Authors: Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
      Purpose: Data structure definition for the problem 
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





#define LINE_BUF_LEN 100

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

#ifndef P_DATA_STRUCTURES
#define P_DATA_STRUCTURES

struct point
  {
    double x;
    double y;
    double prob;
  };

typedef struct
  {
    char*         name;      	 /* instance name */
    char         edge_weight_type[LINE_BUF_LEN];  /* selfexplanatory */
    long int      n;                      /* number of cities */
    struct point  *nodeptr;               /* array of structs containing coordinates of nodes */
    LS_DISTANCE      **distance;	        	/* distance matrix: distance[i][j] gives distance
        													   between city i und j */
  }
problem;


#endif



