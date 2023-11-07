dset ^../data/ERA5_crossspectrum_data.dat
title monthly data regressed to PCs
undef -99999.
xdef 61 linear 0.00 0.00826446281
ydef 1  linear 9.75 0.5
zdef 10 linear 1 1
TDEF 43 linear 00:00Z01JAN1979 1yr
vars 6
re  10 99 real part
im  10 99 imag part
zz  10 99 z power spectrum
mm  10 99 m power spectrum
za  10 99 z autocorrelation
ma  10 99 m autocorrelation
endvars
