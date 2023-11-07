dset ^../data/zm_regression_data.dat
title monthly data regressed to MEI
undef -99999.
xdef 2   linear 0.3125  0.625
ydef 1   linear 9.75 0.5
zdef 4   linear 1 1
TDEF 10  linear 00:00Z01JAN1979 1yr
vars 4
z1 4 99 z1
z2 4 99 z2
m1 4 99 m1
m2 4 99 m2
endvars
