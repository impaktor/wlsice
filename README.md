# Least Squares Including Correlation in Error

Code for performing the weighted Least Squares Including Correlation in Error
(LS-ICE) algorithm when fitting a function to ensemble averages, implemented
in three languages, with example data with M=200 fractional Brownian motion
trajectories (generated with parameters as in paper "Fitting a function to
time-dependent ensemble averaged data"):

1. Python scripts that read in the data and apply the LS-ICE algorithm to it.

2. Octave/Matlab - code ported from Python. `f.m`, `df.m`, and `d2f.m`
   defines the analytical function to fit, its gradient and hessian,
   respectively. Run `example.m` file for example.

3. Common Lisp - Depends on GSLL (GSL Lisp Library), which in turn depends on
   GSL (GNU Scientific library), and libffi. Evaluate example.lisp. Tested on
   SBCL.


# Generating trajectories

We also provide the scripts (Python) used to generate trajectories for the
four example systems.


# License

All code is made available under GNU General Public License v3.
https://www.gnu.org/licenses/gpl-3.0.html
