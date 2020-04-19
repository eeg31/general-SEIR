# general-SEIR
raw generalized SEIR model simulation package with age structure (for Glennon et al. 2019 Phil Trans B)

Provides simulation function for Glennon EE, Becker DJ, Peel AJ, Garnier R, Suu-Ire RD, Gibson L,
Hayman DTS, Wood JLN, Cunningham AA, Plowright RK, and Restif O, "What is stirring in the reservoir? 
Modelling mechanisms of henipavirus circulation in fruit bat hosts," Philosophical Transactions of
the Royal Society B 374(1782), https://doi.org/10.1098/rstb.2019.0021.

Implements Gillespie stochastic simulation algorithm (gSSA) with adaptive tau-leaping to simulate any 
of 46 SEIR-derived compartmental infection models on a population with a 2-parameter seasonal birth pulse, 
three age classes, and waning maternal immunity. Also provides a helper function (getBeta) for converting 
R0 to a transmission rate (as a function of all other model parameters). Parameters are described in the
corresponding paper.
