This software package provides a reasonably high performing implementation of estimation-based local search (iterative improvement) algorithm to tackle the Probabilistic Traveling Salesman Problem (PTSP).


The PTSP is a paradigmatic example of a stochastic combinatorial optimization problem. Estimation-based local search (2.5-opt-EEais) is currently the state-of-the-art iterative improvement algorithm for the PTSP that starts from some initial solution and repeatedly tries to move from a current solution to a lower cost neighboring one. The search terminates in a local optimum, that is, a solution that does not have any improving neighbor. A peculiarity of 2.5-opt-EEais is that the cost of the neighboring solutions are estimated using delta evaluation, a technique that considers only the cost contribution of solution components that are not common between two neighbor solutions. The high performance of this algorithm can be attributed to the adoption of the 2.5-exchange neighborhood relation, neighborhood reduction techniques (fixed-radius search, candidate lists, and don't look bits), and the variance reduction techniques such as the method of common random numbers, adaptive sample size, and importance sampling. For a more detailed explanation, see:

# P. Balaprakash, M. Birattari, T. Stützle, and M. Dorigo. Adaptive sample size and importance sampling in estimation-based local search for the probabilistic traveling salesman problem. European Journal of Operational Research, 199(1):98-110, 2009. [DOI](DOI.md)

# M. Birattari, P. Balaprakash, T. Stützle, and M. Dorigo. Estimation-based local search for stochastic combinatorial optimization using delta evaluations: A case study in the probabilistic traveling salesman problem. INFORMS Journal on Computing, 20(4):644-658, 2008.


Note:
Developed under Linux; Implementation in C; Uses GNU Scientific Library





![http://www.gnu.org/graphics/gplv3-127x51.png](http://www.gnu.org/graphics/gplv3-127x51.png)