dset ^../data/ERA5_z_m_ori_data.dat
title monthly data regressed to PCs
undef -99999.
xdef 1   linear 0.3125  0.625
ydef 1   linear 9.75 0.5
zdef 10  linear 1 1
TDEF 5082  linear 00:00Z01DEC1979 1dy
vars 2
z  10 99 real part
m  10 99 imag part
endvars
