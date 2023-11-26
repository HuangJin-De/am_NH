dset ^../data/ERA5_spectrum_mean.dat
title era5 monthly mean data from 1979 to 2021
undef -99999.
xdef 1 linear 0.3125  0.625
ydef 161 linear 9.75 0.5
zdef 6   levels  1000. 925. 850. 700. 500. 250.
TDEF 365 linear 00:00Z01JAN1979 1dy
vars 4
u 6 99 zonal wind
v 6 99 meridional wind
emc 6 99 eddy momentum flux
mt 1 99 mountain torque
endvars
