/*
#####################################################################
##########    Estimation-based Local Search for the PTSP    #########
#####################################################################

      Version: 1.0
      File:    stopwatch.c
      Authors: Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
      Purpose: Functions for timer  
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
#include <errno.h>
#include <error.h>
#include <sys/time.h>

#define MAX_TIME 36000

void
stopwatch_start(double total_time)
{
  struct itimerval time;

  if ((total_time>MAX_TIME)||(total_time<0))
    error(EXIT_FAILURE,0,"stopwatch: Invalid time");

  total_time = (total_time==0)?MAX_TIME:total_time;

  time.it_value.tv_sec = time.it_interval.tv_sec =
                           (long int)(total_time);
  time.it_value.tv_usec = time.it_interval.tv_usec =
                            (long int)(1e6*(time.it_value.tv_sec-total_time));

  if (setitimer(ITIMER_PROF, &time, 0))
    error(EXIT_FAILURE,errno,"stopwatch: error while setting the timer");
}

double
stopwatch_read()
{
  struct itimerval time;

  if(getitimer(ITIMER_PROF,&time))
    error(EXIT_FAILURE,errno,"stopwatch: error while reading the timer");

  return 0.00001+time.it_interval.tv_sec - time.it_value.tv_sec +
         ((double)(time.it_interval.tv_usec - time.it_value.tv_usec))/1e6;
}

