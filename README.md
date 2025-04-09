# Sensitivity Analysis - Rift Valley Fever Virus

This analysis and code was started by Sophie North (University of Surrey) in 2023, and developed further and completed by Jessica Furber in 2025 for the sensitivity analysis of four parameters in a RVFV model developed by Lo. Iacono et al. [1], when the virus is not present.

The four parameters considered are:

- $\mathcal{A}$: Typical area scanned by mosquitoes flyers looking for oviposition sites.
- $\kappa$: The proportion of the water body area/soil that female adult mosquitoes, lays eggs on.
- $b$: The number of eggs laid per batch by female adults. 
- $\rho^{\max}$: The maximum density of eggs in the water body for mosquitoes to lay eggs.

Further, we firstly consider constant values for the air temperature, water body areas and livestock numbers (values that are brought into the model using empirical data). Though we do vary the values of air temperature, water body areas, and a parameter in the model denoted as $p_f$ (detection probability) for a low, medium and high value to compare the effects on the four sensitive parameters we are investigating. Secondly, we consider periodic functions for the air temperature and water body areas.

The code has been created in R, version 4.4.3.

[1] G. L. Iacono, A. A. Cunningham, B. Bett, D. Grace, D. W. Redding, and J. L. N. Wood, “Environmental limits of Rift Valley fever revealed using ecoepidemiological mechanistic models,” Proceedings of the National Academy of Sciences of the United States of America, vol. 115,pp. E7448–E7456, July 2018.
