dset ^../data/int_data.dat
title monthly data regressed to MEI
undef -99999.
xdef 1   linear 0.3125  0.625
ydef 161 linear 9.75 0.5
zdef 141 linear 1 1
TDEF 42  linear 00:00Z01JAN1979 1yr
vars 3
u    141 99 zonal wind
emc  141 99 eddy
p    141 99 intp
endvars
