# CGPE
This is part of my work on the numerical solution of the Complex Gross-Pitaevskii equation.

You can see the annual presentation.pdf for a overview of this work.

Some code is in the folder CGPE. The CGPEbvp.m use MATLAB function bnp5c to solve the stationary radial solution. And the dynamic_pd.m implement Strang-splitting spectral method to investigate the dynamics of this system. To implement dynamic_pd.m, you should load those data file(xcol,ycol,reycol,imycol) to get the initial solution. You can also get those data by CGPEbvp.m in different parameters.
