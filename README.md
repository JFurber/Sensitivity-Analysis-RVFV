# Sensitivity Analysis - Rift Valley Fever Virus

This analysis and code was started by Sophie North (University of Surrey) in 2023, and developed further and completed by Jessica Furber in 2025 for the sensitivity analysis of a total of ten parameters in a RVFV model developed by Lo. Iacono et al. [1] - four parameters for each of the two mosquito species of interest (\textit{Culex} and \textit{Aedes}) and two parameters shared across species - when the virus is not present.

The main six parameters considered are:

- $\mathcal{A}$: Typical area scanned by mosquitoes flyers looking for oviposition sites.
- $\kappa$: The proportion of the water body area/soil that female adult mosquitoes, lays eggs on.
- $b$: The number of eggs laid per batch by female adults. 
- $\rho$: The maximum density of eggs in the water body for mosquitoes to lay eggs.
- $q$: Parameter for the impact of the livestock on vector fecundity and gonotrophic cycles (or biting rate)
- Livestock Total

The first four parameters are repeated twice for the two different mosquito species considered, with the final two parameters shared across species, as stated.

We firstly consider constant values for the air temperature, water body areas and livestock numbers (values that are brought into the model using empirical data). Though we do vary the values of air temperature, water body areas, and a parameter in the model denoted as $p_f$ (detection probability) for a low, medium and high value to compare the effects on the four sensitive parameters we are investigating (Constant Analysis folder).

Secondly, we consider periodic functions for the air temperature and water body areas, completing a time-varying analysis (Periodic Analysis folder).

The code has been created in R, version 4.4.3.

[1] G. L. Iacono, A. A. Cunningham, B. Bett, D. Grace, D. W. Redding, and J. L. N. Wood, “Environmental limits of Rift Valley fever revealed using ecoepidemiological mechanistic models,” Proceedings of the National Academy of Sciences of the United States of America, vol. 115,pp. E7448–E7456, July 2018.
