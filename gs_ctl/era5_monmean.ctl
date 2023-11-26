dset ^../data/ERA5_spectrum_mean.dat
options yrev
title era5 monthly mean data from 1979 to 2021
undef -99999.
xdef 1    linear 0.3125  0.625
ydef 108  linear 9.75 0.75
zdef 27   levels 1000 975 950 925 900 875 850 825 800 775 750 700 650 600 550 500 450 400 350 300 250 225 200 175 150 125 100
TDEF 365  linear 00:00Z01JAN1979 1mo
vars 4
u   27 99 zonal wind
v   27 99 meridional wind
emc 27 99 eddy momentum flux
mt   1 99 mountain torque
endvars
