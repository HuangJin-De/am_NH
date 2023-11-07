dset ^../data/ERA5_daily_reg_mei.dat
title monthly data regressed to MEI
undef -99999.
xdef 1   linear 0.3125  0.625
ydef 161 linear 9.75 0.5
zdef 6   levels  1000. 925. 850. 700. 500. 250.
TDEF 1 linear 00:00Z01JAN1979 1dy
vars 1
u  6 99 zonal wind
endvars
