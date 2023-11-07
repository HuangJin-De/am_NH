dset ^../data/test.dat
title monthly data regressed to MEI
undef -99999.
xdef 1   linear 0.3125  0.625
ydef 1   linear 9.75 0.5
zdef 1   levels  1000.
TDEF 121 linear 00:00Z01JAN1979 1dy
vars 2
z  1 99 zonal wind
m  1 99 eddy
endvars
