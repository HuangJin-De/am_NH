dset ^../data/ERA5_pc.dat
title monthly data regressed to MEI
undef -99999.
xdef 1   linear 0.3125  0.625
ydef 161 linear 9.75 0.5
zdef 1   levels  1000. 925. 850. 700. 500. 250.
TDEF 10 linear 00:00Z01JAN1979 1dy
vars 1
u  1 99 zonal wind
endvars
