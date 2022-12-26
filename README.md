# STATS534

This final project involves parallel programming in C/C++,R 
that sample parameter beta from a posterior distribution

1) Posterior mode is estimated via Newton Ralphson root-finding algorithm
2) Posterior constant P(D) is estimated by laplace approximation. (INLA is preferred in the case of R)
3) Sampling from the posterior distribution is conducted via MCMC, Metropolis-Hastings Algorithm. (One may use rejection algorithm)
4) Calculate the likelihood by monte-carlo integration
