/*
#####################################################################
##########    Estimation-based Local Search for the PTSP    #########
#####################################################################

      Version: 1.0
      File:    README
      Authors: Prasanna Balaprakash, Mauro Birattari, and Thomas Stuetzle
      Purpose: Instructions for the software usage 
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



This is the README file to the software package ELS-PTSP.

This software package was developed by Prasanna Balaprakash in
connection with his research articles:

# P. Balaprakash, M. Birattari, T. Stuetzle, and M. Dorigo (2008):
  Adaptive sample size and importance sampling in the estimation-based
  local search for the probabilistic traveling salesman
  problem. (European Journal of Operational Research)


# M. Birattari, P. Balaprakash, T. Stuetzle, and M. Dorigo (2008):
  Estimation-based local search for stochastic combinatorial
  optimization using delta evaluations: A case study in the
  probabilistic traveling salesman problem. (INFORMS Journal of
  Computing)
 

The software package is freely available subject to the GNU General Public License v3.

If you use ELS-PTSP in your research, I would appreciate a citation in your publication(s).

This software package provides a reasonably high performing
implementation of estimation-based local search (iterative
improvement) algorithm to tackle the Probabilistic Traveling Salesman
Problem (PTSP).


=========
SUMMARY
=========

The PTSP is a paradigmatic example of a stochastic combinatorial
optimization problem. Estimation-based local search (2.5-opt-EEais) is
currently the state-of-the-art iterative improvement algorithm for the
PTSP that starts from some initial solution and repeatedly tries to
move from a current solution to a lower cost neighboring one. The
search terminates in a local optimum, that is, a solution that does
not have any improving neighbor. A peculiarity of 2.5-opt-EEais is
that the cost of the neighboring solutions are estimated using delta
evaluation, a technique that considers only the cost contribution of
solution components that are not common between two neighbor
solutions. The high performance of this algorithm can be attributed to
the adoption of the 2.5-exchange neighborhood relation, neighborhood
reduction techniques (fixed-radius search, candidate lists, and don't
look bits), and the variance reduction techniques such as the method
of common random numbers, adaptive sample size, and importance
sampling.



=========
CONTENTS
=========


This package consists of the following two directories: 

1) Bin

All the source files of the estimation local search are placed here.

2) Evaluator

All the source files to evaluate the cost of the solution found are placed here.

(More explantion is given in the section USAGE)

Inside Bin/
-----------

The main control routines, including the command line parser:
ptspls.c

Functions to read the PTSP file:
readFile.h
readFile.c

Functions related to the estimation procedure:
sampleLS.h
sampleLS.c

Functions to generate initial solutions based on some heuristics:
heuristics.h
heuristics.c

Functions related to the importance sampling and adaptive sample size procedure:
adaptiveSampling.h
adaptiveSampling.c  

Time measurement:
stopwatch.h
stopwatch.c

The data structure of the problem:
problemdataStructures.h  

The data structure used by the local search:
sampleLSdataStructures.h  

The look up table of the t-test used by the adaptive sample size:
statsTables.h

Some auxillary functions:
sampleLSauxiliaryFunctions.h


Makefile

Note that the instances should be in PTSPLIB format. See the file ch01000-0000001103-0.100.ptsp for an example.
More instances can be obtained from the following URL: http://iridia.ulb.ac.be/supp/IridiaSupp2008-010/


=====
Code
=====


The software was developed in C under Linux, using the GNU 3.3
gcc compiler and extensively tested in this environment. It should be
noted that the software uses GNU/GSL libraries. The software is
distributed as a gzipped tar file.

Go inside Bin and Evaluator directories and type 'make' to compile under Linux; the executables 'ptspls' and 'evaluate' are produced.

======
USAGE
======

Type the following to see the usage:

Bin/ptspls --help

Evaluator/evaluate --help 

Bin/ptspls ch01000-0000001103-0.100.ptsp > els.output.txt 

Evaluator/evaluate ch01000-0000001103-0.100.ptsp els.output.txt > els.results.txt




