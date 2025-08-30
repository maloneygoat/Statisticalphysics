# Statistical Physics
A collection of projects looking at some common models from statistical physics.

# MonteCarlo-2DIsing
This Ising model is one of the most studied models in all of modern physics. The classic model consists of a lattice of specific number of sites, each with a spin value that can be up, +1 or down, -1 and an interaction between them of strength J. The model we will be looking at is the 2D Ising model for a N$\times$N lattice. The Hamiltonian is 

$$H = J\sum_{\langle ij\rangle}s_is_j$$

For our simulations we will take J and $k_B$ to be 1 and the summation is over all nearest neighbours and $s_i, s_j$ are our spins. This model is much more interesting than the simple 1D due to the phase transition that can occurs at a non zero $T_c \approx 2.269$. 

The first exact solution for the 2D Ising was published by Onsagner and displayes a second-order phase transition. The goal is to model the Ising model with Monte Carlo simulations. There are two widely used computational methods for Monte Carlo simulations of the 2D Ising model, the "Single-Flip" and "Cluster" Monte Carlo. 

## Single Flip/Metropolis Algorithm

The Single Flip or Metropolis Algorithm is the simplest model. First we randomly initialise our N$\times$N lattice with different spins. Next we randomly and uniformly choose from the set of all possible states, a configuration that differs from our current one by a single spin flip. This is called the \textit{priori} probability. This probability is constant and is not weighted by energy. A simple example would be a Monte Carlo simulation for throwing a dice, the \textit{priori} probability would be 1 in 6 for any face. After randomly picking a new state we propose changing to that state. This acceptance/rejection criterion is done using the Metropolis criterion, which is a weighting of states based on the energy difference, $\Delta E$, between our current state and the proposed new state. The probability that a state is accepted is 

$$P = \begin{cases}
1 & \text{if } \Delta E \leq 0,\\
e^{-\Delta E /T}  & \text{if } \Delta E > 0
\end{cases}$$

We then proceed to repeat these steps for a set number of iterations, keeping T constant throughout the iterations. This method is simple and incredibly effecting at high temperatures when spins are disordered, however approaching the critical temperature $T_c$, we see a critical slowing in the algorithm. This coupled with only flipping once spin at a time, leads to exploring the configuration space very slowly, requiring a large amount of iterations to get meaningful statistics.

## Cluster/Swendsen-Wang Algorithm

The second model used is the Cluster or Swendsen-Wang Algorithm. We once again initialise a random configuration of our spins to start. Next we choose a random lattice site, called a seed, and begin a breadth first search starting from this seed. We look at the seeds neighbours, if we encounter any neighbours with the same spin we have a probability q of adding them to our cluster. We repeat this process until we have exhausted all available paths. This results in creating a cluster of same spin sites. We propose a next state by flipping all spins in our cluster, the acceptance/rejection criterion is again the Metropolis criterion. The \textit{priori} in this model is the probability q. By choosing q = $1 - e^{2/T}$ we can force that our new state is always accepted. This model immediately fixes the issues we see with the Single Flip method as it allows for large jumps in configuration space between two steps.

## Simulations
Simulations were performed on 10 $\times$ 10, 20 $\times$ 20 and 100 $\times$ 100 lattices with 10,000 iterations at each temperature to ensure good convergence. As we increase lattice size we can see that our simulation is converging towards the theoretical model. This is to be expected as the Onsarger model assumes infinite lattice size so as we increase our lattice size we expect to see convergence.

# Numerical Modelling of Real Gas
The van der Waals equation describes the he behavior of real gases. It is an equation of state that relates the pressure, volume, number of molecules, and temperature in a fluid. We are going to study the equation and the phase diagram with respect to Mercury. The van der Waals equation is 

$$p = \frac{RT}{V/n - b} - \frac{a}{(V/n)^2}$$ 

where R is the gas constant, n is number of moles which we set to 1. For Mercury a = 8.2 $L^2bar/mol^2$ and b = 0.01696 $L/mol$.

We can also write this equation in the more universal reduced van der Waals:

$$p_r = \frac{8T_r}{3V_r-1} - \frac{3}{V_r^2}$$ 

where $V_r = V/V_c, \text{   } p_r = p/p_c\text{ and }T_r = T/T_c $. $V_c, \text{   } p_c \text{ and }T_c$ are our critical values. These values depend on the constants a, b and are 

$$V_c = 3b, \text{   } p_c = \frac{a}{27b^2}, \text{   } T_c = \frac{8a}{27Rb} $$

## Bindoal Curve: Maxwells Construction
With the above formulae we can plot the isotherms and we can now try to find the binodal line using the equal area method or Maxwells Construction. The binodal curve, also known as the coexistence curve, represents the boundary on a phase diagram where two distinct phases can coexist in equilibrium By drawing a line though our Isotherms for any temperature below $T_c$, if the area enclosed by this above and below this line equal, then the first and last intersection points are on the binodal curve. By repeating this construction for different $T< T_c$ we can acquire enough points to draw this curve. 

This was done numerically by first guessing a line that will intersect a specific isotherm three times. If the line did not intersect 3 times, or the area was not equal we adjusted the line up or down using a root finder algorithm. This process was repeated for each value of $T$ below $T_c$.

## Spinodal Curve
The spinodal curve is much easier to find than the Binodal, requring only simple derivatives. By taking the first derivative of the van der Waals equation we obtain

$$ \left( \frac{\partial p}{\partial V} \right)_T = \frac{a}{V^2} - \frac{RT}{(V-b)^2}$$ 

By setting this equation to 0, solving for T and subbing into our oringinal van der Waals we get 

$$p = \frac{a}{V^2} - \frac{2ab}{V^3}$$ 

which we can once again write in reduced form using our critical values. After subbing in we get the equation

$$ p_r = \frac{3}{V_r^2} - \frac{2}{V_r^3}$$

# 2D random walk on a Square Lattice
Plotting and calculating the mean square end-end distance for a random walk and a self-avoiding random walk. 

 $$\langle R^2\rangle = \frac{1}{\Omega}\sum_{i=1}^{\Omega}R_i^2$$


