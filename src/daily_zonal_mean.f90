program daily_zonal_mean
use netcdf
use, intrinsic :: iso_c_binding
implicit none

include 'fftw3.f03'

integer, parameter :: nx=480,ny=241,nz=37
real, parameter :: tri_pi=4.*atan(1.)
real, dimension(nz) :: lev
real, dimension(ny) :: lat
real, dimension(nx) :: lon
integer :: i,j,k,m,n,o,t
integer :: nt,sny,snz,tnt
integer :: ierr,ncid1,varid1
integer :: yr,da,access
integer :: yrs,yre
integer, dimension(4) :: start4,count4
integer, dimension(3) :: start3,count3
integer, dimension(12) :: mo_da
real :: lats,late
real :: dum1,dum2,dum3,dum4,dum5
real, dimension(:,:,:,:), allocatable :: u,v,emc
real, dimension(:,:,:), allocatable :: um,vm,emcm
real, dimension(:,:,:), allocatable :: mt
real, dimension(:,:), allocatable :: mtm,mta, ua,va
character(200) :: path,fname,tdum1
real(c_double), dimension(:), allocatable :: fin
complex(c_double_complex), dimension(:), allocatable :: fcoe
TYPE(c_ptr) :: plan_forward, plan_backward

path="/data/der0318/work/am_NH/"

lats=10
late=90
yrs=1979
yre=2017

yr=yrs
write(fname,'(2A,I4,A)') trim(path),"/ERA-I/U/daily_interim_U_",yr,".nc"

ierr=nf90_open(trim(fname),nf90_nowrite,ncid1)
if (ierr/=nf90_noerr) write(*,*) "open fail"
ierr=nf90_inq_varid(ncid1,'latitude',varid1)
if (ierr/=nf90_noerr) write(*,*) "inq var fail"
ierr=nf90_get_var(ncid1,varid1,lat)
if (ierr/=nf90_noerr) write(*,*) "read fail"
ierr=nf90_inq_varid(ncid1,'level',varid1)
if (ierr/=nf90_noerr) write(*,*) "inq var fail"
ierr=nf90_get_var(ncid1,varid1,lev)
if (ierr/=nf90_noerr) write(*,*) "read fail"
ierr=nf90_close(ncid1)
if (ierr/=nf90_noerr) write(*,*) "close fail"

lev=lev*100.

tnt=0
do yr=yrs,yre
  tnt=tnt+365
  if (mod(yr,4)==0) tnt=tnt+1
enddo
write(*,*) tnt

i=minloc(abs(lat-lats),1)
j=minloc(abs(lat-late),1)
if (i>=j) then
  k=i
  i=j
  j=k
endif
sny=j-i+1

write(*,*) lat(i),lat(j),sny

start4=(/1,i,1,1/)
start3=(/1,i,1/)

snz=minloc(abs(lev-10000.),1)

allocate(um(sny,snz,365),vm(sny,snz,365),emcm(sny,snz,365),mtm(sny,365))
allocate(ua(sny,snz),va(sny,snz))
um=0.; vm=0.; emcm=0.; mtm=0.;

n=1
do yr=yrs,yre
  mo_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  if (mod(yr,4)==0) mo_da(2)=mo_da(2)+1
  nt=sum(mo_da,1)
  count4=(/nx,sny,snz,nt/)
  count3=(/nx,sny,nt/)
  allocate(u(nx,sny,snz,nt),v(nx,sny,snz,nt),emc(nx,sny,snz,nt),mt(nx,sny,nt))

  write(fname,'(2A,I4,A)') trim(path),"/ERA-I/U/daily_interim_U_",yr,".nc"
  !write(*,*) trim(fname),access(trim(fname),' ')
 
  ierr=nf90_open(trim(fname),nf90_nowrite,ncid1)
  if (ierr/=nf90_noerr) write(*,*) "open fail"
  ierr=nf90_inq_varid(ncid1,'u',varid1)
  if (ierr/=nf90_noerr) write(*,*) "inq var fail"
  ierr=nf90_get_var(ncid1,varid1,u,start=start4,count=count4)
  if (ierr/=nf90_noerr) write(*,*) "read fail"
  ierr=nf90_close(ncid1)
  if (ierr/=nf90_noerr) write(*,*) "close fail"
  write(*,*) "read u"

  write(fname,'(2A,I4,A)') trim(path),"/ERA-I/V/daily_interim_V_",yr,".nc"
  !write(*,*) trim(fname),access(trim(fname),' ')

  ierr=nf90_open(trim(fname),nf90_nowrite,ncid1)
  if (ierr/=nf90_noerr) write(*,*) "open fail"
  ierr=nf90_inq_varid(ncid1,'v',varid1)
  if (ierr/=nf90_noerr) write(*,*) "inq var fail"
  ierr=nf90_get_var(ncid1,varid1,v,start=start4,count=count4)
  if (ierr/=nf90_noerr) write(*,*) "read fail"
  ierr=nf90_close(ncid1)
  if (ierr/=nf90_noerr) write(*,*) "close fail"
  write(*,*) "read v"

  write(fname,'(2A,I4,A)') trim(path),"/ITM/MT/MT_",yr,".nc"

  ierr=nf90_open(trim(fname),nf90_nowrite,ncid1)
  if (ierr/=nf90_noerr) write(*,*) "open fail"
  ierr=nf90_inq_varid(ncid1,'mt',varid1)
  if (ierr/=nf90_noerr) write(*,*) "inq var fail"
  ierr=nf90_get_var(ncid1,varid1,mt,start=start3,count=count3)
  if (ierr/=nf90_noerr) write(*,*) "read fail"
  ierr=nf90_close(ncid1)
  if (ierr/=nf90_noerr) write(*,*) "close fail"
  write(*,*) "read mt"!,maxval(mt)


  do t=1,nt
    o=t
    if (mod(yr,4)==0 .and. t==60) goto 444
    if (mod(yr,4)==0 .and. t>60) o=t-1

    ua(:,:)=sum(u(:,:,:,t),1)/real(nx)
    va(:,:)=sum(v(:,:,:,t),1)/real(nx)
    do i=1,nx
      emc(i,:,:,t)=(u(i,:,:,t)-ua(:,:))*(v(i,:,:,t)-va(:,:))
    enddo

    um(:,:,o)=um(:,:,o)+sum(u(:,:,:,t),1)/real(nx)/real(yre-yrs+1)
    vm(:,:,o)=vm(:,:,o)+sum(v(:,:,:,t),1)/real(nx)/real(yre-yrs+1)
    emcm(:,:,o)=emcm(:,:,o)+sum(emc(:,:,:,t),1)/real(nx)/real(yre-yrs+1)
    mtm(:,o)=mtm(:,o)+sum(mt(:,:,t),1)/real(nx)/real(yre-yrs+1)
     
    444 continue
  enddo

  deallocate(u,v,emc,mt)
  n=n+nt
  write(*,*) yr
enddo

write(*,*) "after mean"

! calculate mean using fft
i=365
allocate(fin(i),fcoe(i/2+1))
plan_forward=fftw_plan_dft_r2c_1d(i,fin,fcoe,FFTW_ESTIMATE)
plan_backward=fftw_plan_dft_c2r_1d(i,fcoe,fin,FFTW_ESTIMATE)

do k=1,snz
do j=1,sny
  fin=um(j,k,:)
  call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
  fcoe(5:i/2+1)=0.
  call fftw_execute_dft_c2r(plan_backward,fcoe,fin)
  um(j,k,:)=fin/real(i)

  fin=vm(j,k,:)
  call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
  fcoe(5:i/2+1)=0.
  call fftw_execute_dft_c2r(plan_backward,fcoe,fin)
  vm(j,k,:)=fin/real(i)

  fin=emcm(j,k,:)
  call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
  fcoe(5:i/2+1)=0.
  call fftw_execute_dft_c2r(plan_backward,fcoe,fin)
  emcm(j,k,:)=fin/real(i)

  if (k==1) then
    fin=mtm(j,:)
    call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
    fcoe(5:i/2+1)=0.
    call fftw_execute_dft_c2r(plan_backward,fcoe,fin)
    mtm(j,:)=fin/real(i)
  endif
enddo
enddo

call fftw_free(plan_forward)
call fftw_free(plan_backward)


write(*,*) "after fft mean"!,mtm

fname=trim(path)//"/data/ERA5_spectrum_mean.dat"
open(10,file=trim(fname),access="direct",recl=sny*snz*3+sny)
do t=1,365
  write(10,rec=t) um(:,:,t),vm(:,:,t),emcm(:,:,t),mtm(:,t)
enddo
close(10)

end program daily_zonal_mean
