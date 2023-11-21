# compiler
FC = mpifort -fpp -DPURE
DEBUG = -CB -g -traceback -check all,noarg_temp_created -debug all
FCFLAGS = -O3 -free -mcmodel=large -heap-arrays 10 -shared-intel -fp-model precise #$(DEBUG)
FINCLUDE = -I/home/der0318/.local/include
LDLIBS = -L/home/der0318/.local/lib -lfftw3 -lnetcdff -lnetcdf -lhdf5 -lhdf5_hl -lsz -Wl,-rpath,/home/der0318/.local/lib -L/home/der0318/.local/lib -lflapack -lfblas

# code paths
VPATH = src

# objects
#LIST = month_mean.f90 
#a.out: month_mean.o
#month_mean.o: month_mean.f90

#LIST = daily_zonal_mean.f90
#a.out: daily_zonal_mean.o
#daily_zonal_mean.o: daily_zonal_mean.f90

#LIST = remove_enso_monthly.f90
#a.out: remove_enso_monthly.o
#remove_enso_monthly.o: remove_enso_monthly.f90

#LIST = remove_enso_daily.f90
#a.out: remove_enso_daily.o
#remove_enso_daily.o: remove_enso_daily.f90

#LIST = projection_onto_eof.f90
#a.out: projection_onto_eof.o
#projection_onto_eof.o: projection_onto_eof.f90

#LIST = spectrum_zm.f90
#a.out: spectrum_zm.o
#spectrum_zm.o: spectrum_zm.f90

LIST = retopy.f90
a.out: retopy.o
retopy.o: retopy.f90

#LIST = calculate_mountain_torque.f90
#a.out: calculate_mountain_torque.o
#calculate_mountain_torque.o: calculate_mountain_torque.f90

#LIST =  pca_analysis.f90 forpca.f90 forsvd.f90 foreig.f90 kinds.f90
#a.out: pca_analysis.o
#pca_analysis.o: pca_analysis.f90 forpca.o kinds.o
#forpca.o: forpca.f90 kinds.o foreig.o forsvd.o
#forsvd.o: forsvd.f90 kinds.o
#foreig.o: foreig.f90 kinds.o
#kinds.o: kinds.f90

LIST_o = $(LIST:.f90=.o)
target = a.out

all: $(target)

$(LIST_o): %.o: %.f90
	$(FC) $(FCFLAGS) $(FINCLUDE) -c $< 

$(target) : $(LIST_o)
	$(FC) $(FCFLAGS) $(FINCLUDE) $^ -o $@ $(LDLIBS)

clean:
	rm -rf *.o *.mod a.out
