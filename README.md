Supplemental document:
## Exploration of geologically reasonable surface uplift and horizontal propagation rates
Kleber, et al.

# Overview

Reconstruction of displaced landforms with their ages provides the basic information to determine rates of faulting or fold growth. Gold and Cowgill (2011) and numerous others have explored this issue by considering the age-displacement pairs as part of a range of possible histories constrained by plausible samples through the age and position ranges. This approach has advantages over simple linear fits through the data including assessment of secular variations in rates, exploration of the age and position constraints, and the general nuance of the implications of the resulting history on the processes under investigation.

In the case of Wheeler Ridge, the deformed surfaces provide information about surface uplift rate and horizontal propagation rates. We assume that the ages are normally distributed about a mean by a specified standard deviation (). We truncate the input range at 2 We also assume that the surface uplift (SU) and horizontal position (propagation—displacement zone of the fold tip) are uniformly distributed between the specified ranges of preservation of the surfaces along the Wheeler Ridge axis. The age-SU and age-propagation pairs are selected from the specified ranges by random (“Monte Carlo”) sampling (n = 105-106 times). We subject to a positivity constraint, and include only those histories with positive rates of uplift or propagation. The resulting sequential histories from surface to surface are concatenated to develop a cumulative rate distribution which is then trimmed for significance. The analysis is done in Matlab and the code is available here.
