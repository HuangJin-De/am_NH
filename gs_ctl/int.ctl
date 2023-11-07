dset ^../data/int_data.dat
title monthly data regressed to MEI
undef -99999.
xdef 1   linear 0.3125  0.625
ydef 161 linear 9.75 0.5
zdef 121 linear 1 1
TDEF 42  linear 00:00Z01JAN1979 1yr
vars 3
u    121 99 zonal wind
emc  121 99 eddy
p    121 99 intp
endvars
