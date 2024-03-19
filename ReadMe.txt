Hello!

This is a brief overview of the tyre fitting code, along with assumptions made.

The PAC2002 equations have been used, based on the ADAMS tyre help file. MF5.2 was not used due to the reduced dependency of longitudinal forces on camber for the analysed tyres. 

FMINCON optimiser has been used, with non-linear constraints to ensure that resultant coefficients are within the range of values stipulated by literature. These values and relations for each case (pure lateral, longitudinal) have been stated in the functions that have the prefix 'nonlcon', as can be found within the folder (MF52 Equations).

The steps to run are as follows:
1. Please run the function named MF52_MASTER_SCRIPT.m
2. Select lateral data file first, I have included one called 'lateral_data.m'
3. Lateral fitting process will be completed with the coefficient results automatically saved
4. Select longitudinal data file, I have included one called 'longitudinal_data.m'
5. Here onwards, the code will run on its own to the end, at which point you will receive the individual and combined slip modelling results

Known errors:
1. FMINCON solver has not been optimised thoroughly, as a result, there are warnings of matrices working with singular precision
2. Combined slip coefficient constraints led to instability and infeasible results at time, hence they have not been included.
3. Inconsistency in LFZ0 value for combined model, 1.1 was selected for fitting process (PAC2002_FY_Combined.m & PAC2002_FX_Combined.m) and 1.2 is used for output in 'MF52_Combined.m'. This is due to the maximum nominal load tested being below the expected wheel loads for OBR23, and use LFZ0 of 1 resulted in poor fitting as expected when vertical loads were extrapolated. 

Thank you for reading!
- Ananth