Solomon Polachek and Bong Joon Yoon, "Panel Estimates of a Two-Tiered
Earnings Frontier", Journal of Applied Econometrics, Vol. 11, No. 2, 
1996, pp. 169-178.

Data Source:  Panel Study of Income Dynamics 1969-84

The data set is an ASCII file, PSIDQ.DAT (zipped in py-data.zip).  The
data are in Fortran format:

YEAR WAGE EDUC EXP TEN SEQN68
I2, F10.5, 4I5

The variables used in the study are:
YEAR (69 - 84), log of WAGE in 1967 constant Dollars,
EDUCation and EXPerience in years, TENure in months;
SEQN68 is person identification number.

Number of Observations = 13408 (lines) = 838 ind x 16 years
