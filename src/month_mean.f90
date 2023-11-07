program month_mean
use netcdf
use, intrinsic :: iso_c_binding
implicit none

include 'fftw3.f03'

integer, parameter :: nx=576,ny=360,nz=6
real, parameter :: tri_pi=4.*atan(1.)
real, dimension(nz), parameter :: lev=(/1000.,925.,850.,700.,500.,250./)
real, dimension(ny) :: lat
real, dimension(nx) :: lon
integer :: i,j,k,m,n,o,t
integer :: nt,sny,tnt,mnt
integer :: ierr,ncid1,varid1
integer :: yrs,yre
integer :: yr,da,access
integer, dimension(3) :: start,count
integer, dimension(12) :: mo_da
real :: lats,late
real :: dum1,dum2,dum3,dum4,dum5
real, dimension(:,:,:,:), allocatable :: u
real, dimension(:,:,:), allocatable :: um, uc
character(200) :: path,fname,tdum1
real(c_double), dimension(:), allocatable :: fin
complex(c_double_complex), dimension(:), allocatable :: fcoe
TYPE(c_ptr) :: plan_forward, plan_backward

path="/work/der0318/work/am_NH/"

lats=10
late=90
yrs=1979
yre=2021

tnt=(yre-yrs+1)*12

yr=yrs
k=1
write(tdum1,'(I10)') int(lev(k))
write(fname,'(4A,I4,A)') trim(path),"/ERA5/u/u",trim(adjustl(tdum1)),"_",yr,".nc"

ierr=nf90_open(trim(fname),nf90_nowrite,ncid1)
if (ierr/=nf90_noerr) write(*,*) "open fail"
ierr=nf90_inq_varid(ncid1,'lat',varid1)
if (ierr/=nf90_noerr) write(*,*) "inq var fail"
ierr=nf90_get_var(ncid1,varid1,lat)
if (ierr/=nf90_noerr) write(*,*) "read fail"
ierr=nf90_close(ncid1)
if (ierr/=nf90_noerr) write(*,*) "close fail"

i=minloc(abs(lat-lats),1)
j=minloc(abs(lat-late),1)
sny=j-i+1

write(*,*) lat(i),lat(j),sny

start=(/1,i,1/)

allocate(um(sny,nz,tnt),uc(sny,nz,12))

n=1
do yr=yrs,yre
  mo_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  if (mod(yr,4)==0) mo_da(2)=mo_da(2)+1
  nt=sum(mo_da,1)
  count=(/nx,sny,nt/)
  allocate(u(nx,sny,nz,nt))

  do k=1,6
    write(tdum1,'(I10)') int(lev(k))
    write(fname,'(4A,I4,A)') trim(path),"/ERA5/u/u",trim(adjustl(tdum1)),"_",yr,".nc"
    !write(*,*) trim(fname),access(trim(fname),' ')
 
    ierr=nf90_open(trim(fname),nf90_nowrite,ncid1)
    if (ierr/=nf90_noerr) write(*,*) "open fail"
    ierr=nf90_inq_varid(ncid1,'u',varid1)
    if (ierr/=nf90_noerr) write(*,*) "inq var fail"
    ierr=nf90_get_var(ncid1,varid1,u(:,:,k,:),start=start,count=count)
    if (ierr/=nf90_noerr) write(*,*) "read fail"
    ierr=nf90_close(ncid1)
    if (ierr/=nf90_noerr) write(*,*) "close fail"
  enddo

  do i=1,12
    j=sum(mo_da(1:i))-mo_da(i)+1
    k=sum(mo_da(1:i))
    o=n+i-1
    um(:,:,o)=sum(sum(u(:,:,:,j:k),4),1)/real((k-j+1)*nx)
    uc(:,:,i)=uc(:,:,i)+um(:,:,o)/real(yre-yrs+1)
  enddo 

  n=n+12
  deallocate(u)
  write(*,*) yr
enddo
n=n-1
write(*,*) n

! calculate mean using fft
i=12
allocate(fin(i),fcoe(i/2+1))
plan_forward=fftw_plan_dft_r2c_1d(i,fin,fcoe,FFTW_ESTIMATE)
plan_backward=fftw_plan_dft_c2r_1d(i,fcoe,fin,FFTW_ESTIMATE)

do k=1,nz
do j=1,sny
  fin=uc(j,k,:)
  call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
  fcoe(5:i/2+1)=0.
  call fftw_execute_dft_c2r(plan_backward,fcoe,fin)
  uc(j,k,:)=fin/real(i)
enddo
enddo

call fftw_free(plan_forward)
call fftw_free(plan_backward)

deallocate(fin,fcoe)

fname=trim(path)//"/data/ERA5_month_mean_spectrum.dat"
open(10,file=trim(fname),access="direct",recl=sny*nz)
do t=1,12
  write(10,rec=t) uc(:,:,t)
enddo
close(10)

deallocate(um,uc)

end program month_mean
