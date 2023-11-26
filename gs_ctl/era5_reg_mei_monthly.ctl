dset ^../data/ERA5_month_reg_mei.dat
options yrev
title monthly data regressed to MEI
undef -99999.
xdef 1   linear 0.3125  0.625
ydef 108 linear 9.75 0.75
zdef 27   levels 1000 975 950 925 900 875 850 825 800 775 750 700 650 600 550 500 450 400 350 300 250 225 200 175 150 125 100
TDEF 1 linear 00:00Z01JAN1979 1dy
vars 3
u   27 99 zonal wind
emc 27 99 eddy momentum flux
mt   1 99 mountain torque
endvars
