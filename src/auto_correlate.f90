program auto_correlate
use netcdf
implicit none

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
integer :: yrs,yre,ts,te,is,ie
integer, dimension(3) :: start,count
integer, dimension(12) :: mo_da
real :: correlation
real :: lats,late
real :: dum1,dum2,dum3,dum4,dum5,dum6
real, dimension(:), allocatable :: zm,mm,zv,mv,zk,mk
real, dimension(:), allocatable :: hann,x,y
real, dimension(:,:), allocatable :: zp,mp
real, dimension(:,:,:), allocatable :: tmp
character(200) :: path,fname,tdum1

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

allocate(tmp(nt,34,42),x(nt),y(nt))

do t=1,42
  ts=(t-1)*nt+1
  te=t*nt

  do k=1,10
    do i=-60,60
      is=ts+i
      ie=te+i
      if (is<ts) is=ts
      if (ie>te) ie=te
      !if (t/=0 .and. k==1) write(*,*) i, is-ts+1,ie-ts+1,is-i-ts+1,ie-i-ts+1
 
      m=ie-is+1
      x=zp(is:ie,k)
      y=zp(is-i:ie-i,k)
      tmp(61+i,k,t)=correlation(m,x,y)
      if (t/=0 .and. k==1) write(*,*) i,tmp(61+i,k,t)

      x=mp(is:ie,k)
      y=mp(is-i:ie-i,k)
      tmp(61+i,10+k,t)=correlation(m,x,y)

      x=zp(is:ie,k)
      y=mp(is-i:ie-i,k)
      tmp(61+i,20+k,t)=correlation(m,x,y) 
    enddo
  enddo
 
  do i=-60,60
    is=ts+i
    ie=te+i
    if (is<ts) is=ts
    if (ie>te) ie=te

     x=zp(is:ie,1)
     y=mp(is-i:ie-i,2)
     tmp(61+i,31,t)=correlation(m,x,y)
  
     x=zp(is:ie,2)
     y=mp(is-i:ie-i,1)
     tmp(61+i,32,t)=correlation(m,x,y)
  
     x=zp(is:ie,1)
     y=zp(is-i:ie-i,2)
     tmp(61+i,33,t)=correlation(m,x,y)
  enddo

  write(*,*) t
enddo

open(10,file="./data/auto_corr_data.dat",access="direct",recl=nt*34*42)
write(10,rec=1) tmp 
close(10)

end program auto_correlate

real function correlation(n,x,y)

integer, intent(in) :: n
real, dimension(n), intent(in) :: x,y
real :: dum1,dum2,dum3,dum4,dum5

dum1=sum(x,1)/real(n)
dum2=sum(y,1)/real(n)

dum3=sum((x-dum1)*(y-dum2),1)
dum4=sqrt(sum((x-dum1)**2,1))
dum5=sqrt(sum((y-dum2)**2,1))

correlation=dum3/dum4/dum5

end function correlation
