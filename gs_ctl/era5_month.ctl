dset ^../data/ERA5_month_mean.dat
title era5 monthly mean data from 1979 to 2021
undef -99999.
xdef 1 linear 0.3125  0.625
ydef 108 linear 9.75 0.75
zdef 27   levels 1000 975 950 925 900 875 850 825 800 775 750 700 650 600 550 500 450 400 350 300 250 225 200 175 150 125 100
TDEF 172 linear 00:00Z01JAN1979 1mo
vars 1
u 27 99 zonal wind
endvars
