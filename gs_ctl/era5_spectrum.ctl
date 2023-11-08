dset ^../data/ERA5_spectrum_data.dat
options yrev
title era5 spectrum
undef -99999.
xdef 1   linear 0.3125  0.625
ydef 121 linear -60 1
zdef 10  linear 0 1
TDEF 42 linear 00:00Z01JAN1979 1yr
vars 7
z    10 99 z
m    10 99 m
zm   10 99 zm
z1m2  1 99 z1m2
z2m1  1 99 z2m1
z1z2  1 99 z1z2
tmp2  1 99 z2m1
endvars
