program spectrum_zm
use netcdf
use, intrinsic :: iso_c_binding
implicit none

include 'fftw3.f03'

integer, parameter :: nx=576,ny=360,nz=6
real, parameter :: tri_pi=4.*atan(1.),grav=9.8,a=6378000.,d2r=tri_pi/180.
real, dimension(nz), parameter :: lev=(/1000.,925.,850.,700.,500.,250./)
real, dimension(ny) :: lat
real, dimension(:), allocatable :: slat
real, dimension(nx) :: lon
integer :: i,j,k,m,n,o,t
integer :: nt,sny,tnt
integer :: ierr,ncid1,varid1
integer :: yr,da,access
integer :: yrs,yre,ts,te
integer, dimension(3) :: start,count
integer, dimension(12) :: mo_da
real :: lats,late
real :: dum1,dum2,dum3,dum4,dum5,dum6
real, dimension(:), allocatable :: zm,mm,zv,mv,zk,mk
real, dimension(:), allocatable :: hann
real, dimension(:,:), allocatable :: zp,mp,tmp
real, dimension(:,:,:), allocatable :: ab 
character(200) :: path,fname,tdum1

real(c_double), dimension(:), allocatable :: fin,re,im
complex(c_double_complex), dimension(:), allocatable :: fcoe,zn,mn
TYPE(c_ptr) :: plan_forward, plan_backward

path="/work/der0318/work/am_NH/"

lats=10
late=90
yrs=1979
yre=2021

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

tnt=(yre-yrs+1)*4

i=minloc(abs(lat-lats),1)
j=minloc(abs(lat-late),1)
sny=j-i+1

allocate(slat(sny))

slat=lat(i:j)

tnt=0
do yr=yrs,yre
  mo_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  !if (mod(yr,4)==0) mo_da(2)=mo_da(2)+1
  tnt=tnt+sum(mo_da(1:3),1)+mo_da(12)
  if (yr==yrs) tnt=tnt-sum(mo_da(1:3),1)
  if (yr==yre) tnt=tnt-mo_da(12)
enddo
write(*,*) tnt

allocate(zp(tnt,10),mp(tnt,10))

fname="./data/ERA5_z_m_data.dat"
open(10,file=trim(fname),access="direct",recl=20)
do t=1,tnt
  read(10,rec=t) zp(t,:),mp(t,:)
enddo
close(10)

nt=121

allocate(fin(nt),fcoe(nt/2+1),zk(nt),mk(nt),zn(nt/2+1),mn(nt/2+1),hann(nt)) 
allocate(re(nt/2+1),im(nt/2+1))
allocate(zm(1),mm(1),zv(1),mv(1))

dum1=real(nt+1)/2.
do i=1,nt
  !hann(i)=cos(tri_pi*(real(i)-dum1)/real(nt-1))**2 
  hann(i)=0.5-0.5*cos(2*tri_pi*(real(i-1))/real(nt))
  !write(*,*) i, hann(i), cos(tri_pi*(real(i)-dum1)/real(nt))**2
enddo


allocate(tmp(nt,34))

plan_forward=fftw_plan_dft_r2c_1d(nt,fin,fcoe,FFTW_ESTIMATE)
plan_backward=fftw_plan_dft_c2r_1d(nt,fcoe,fin,FFTW_ESTIMATE)


fname="./data/ERA5_spectrum_data.dat"
open(10,file=trim(fname),access="direct",recl=nt*34)

fname="./train_data/ERA5_spectrum_data.dat"
open(20,file=trim(fname),access="direct",recl=nt*34)

n=1

do t=1,42
  ts=(t-1)*nt+1
  te=t*nt

  tmp=0.

  do m=1,10
    zk=zp(ts:te,m)*hann
    mk=mp(ts:te,m)*hann

    zm(1)=sum(zk,1)/real(nt)
    mm(1)=sum(mk,1)/real(nt)
    zv=sqrt(sum((zk-zm(1))**2.,1)/real(nt))
    mv=sqrt(sum((mk-mm(1))**2.,1)/real(nt))

    zk=(zk-zm(1))/zv(1)
    mk=(mk-mm(1))/mv(1)
 
    fin=zk
    call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
    zn=fcoe

    fin=mk
    call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
    mn=fcoe

    fcoe=abs(zn)**2
    call fftw_execute_dft_c2r(plan_backward,fcoe,fin)
    fin=fin/real(nt)**2
  
    tmp(1:nt/2,m)=fin(nt/2+2:nt)
    tmp(nt/2+1:nt,m)=fin(1:nt/2+1)

    fcoe=abs(mn)**2
    call fftw_execute_dft_c2r(plan_backward,fcoe,fin)
    fin=fin/real(nt)**2

    tmp(1:nt/2,10+m)=fin(nt/2+2:nt)
    tmp(nt/2+1:nt,10+m)=fin(1:nt/2+1)

    fcoe=zn*conjg(mn)
    call fftw_execute_dft_c2r(plan_backward,fcoe,fin)
    fin=fin/real(nt)**2

    tmp(1:nt/2,20+m)=fin(nt/2+2:nt)
    tmp(nt/2+1:nt,20+m)=fin(1:nt/2+1)
  enddo

  zk=zp(ts:te,1)*hann
  mk=mp(ts:te,2)*hann

  zm(1)=sum(zk,1)/real(nt)
  mm(1)=sum(mk,1)/real(nt)
  zv=sqrt(sum((zk-zm(1))**2.,1)/real(nt))
  mv=sqrt(sum((mk-mm(1))**2.,1)/real(nt))

  zk=(zk-zm(1))/zv(1)
  mk=(mk-mm(1))/mv(1)

  fin=zk
  call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
  zn=fcoe

  fin=mk
  call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
  mn=fcoe

  fcoe=mn*conjg(zn)
  call fftw_execute_dft_c2r(plan_backward,fcoe,fin)
  fin=fin/real(nt)**2

  tmp(1:nt/2,31)=fin(nt/2+2:nt)
  tmp(nt/2+1:nt,31)=fin(1:nt/2+1)

  zk=zp(ts:te,2)*hann
  mk=mp(ts:te,1)*hann

  zm(1)=sum(zk,1)/real(nt)
  mm(1)=sum(mk,1)/real(nt)
  zv=sqrt(sum((zk-zm(1))**2.,1)/real(nt))
  mv=sqrt(sum((mk-mm(1))**2.,1)/real(nt))

  zk=(zk-zm(1))/zv(1)
  mk=(mk-mm(1))/mv(1)

  fin=zk
  call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
  zn=fcoe

  fin=mk
  call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
  mn=fcoe

  fcoe=mn*conjg(zn)
  call fftw_execute_dft_c2r(plan_backward,fcoe,fin)
  fin=fin/real(nt)**2

  tmp(1:nt/2,32)=fin(nt/2+2:nt)
  tmp(nt/2+1:nt,32)=fin(1:nt/2+1)

  zk=zp(ts:te,2)*hann
  mk=zp(ts:te,1)*hann

  zm(1)=sum(zk,1)/real(nt)
  mm(1)=sum(mk,1)/real(nt)
  zv=sqrt(sum((zk-zm(1))**2.,1)/real(nt))
  mv=sqrt(sum((mk-mm(1))**2.,1)/real(nt))

  zk=(zk-zm(1))/zv(1)
  mk=(mk-mm(1))/mv(1)

  fin=zk
  call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
  zn=fcoe

  fin=mk
  call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
  mn=fcoe

  fcoe=mn*conjg(zn)
  call fftw_execute_dft_c2r(plan_backward,fcoe,fin)
  fin=fin/real(nt)**2

  tmp(1:nt/2,33)=fin(nt/2+2:nt)
  tmp(nt/2+1:nt,33)=fin(1:nt/2+1)

  !write(*,*) t,maxval(fin,1),maxloc(fin,1)
  write(10,rec=n) tmp
  write(20,rec=n) tmp
  n=n+1
enddo


call fftw_free(plan_forward)
call fftw_free(plan_backward)

close(10)
close(20)

deallocate(fin,fcoe)

end program spectrum_zm
