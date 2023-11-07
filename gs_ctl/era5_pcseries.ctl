dset ^../data/ERA5_pc_series.dat
title monthly data regressed to MEI
undef -99999.
xdef 172  linear 0.3125  0.625
ydef 1   linear 9.75 0.5
zdef 1   linear 1 1
TDEF 172 linear 00:00Z01JAN1979 1dy
vars 1
u  1 99 zonal wind
endvars
