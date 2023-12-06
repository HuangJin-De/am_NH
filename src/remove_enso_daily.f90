program daily_zonal_mean
use netcdf
implicit none

integer, parameter :: nx=576,ny=360,nz=6
real, parameter :: tri_pi=4.*atan(1.)
real, dimension(nz), parameter :: lev=(/1000.,925.,850.,700.,500.,250./)
real, dimension(ny) :: lat
real, dimension(nx) :: lon
integer :: i,j,k,m,n,o,t,ii,jj
integer :: nt,sny,tnt
integer :: ierr,ncid1,varid1
integer :: yr,da,access
integer :: yrs,yre
integer, dimension(3) :: start,count
integer, dimension(12) :: mo_da
real :: lats,late
real :: dum1,dum2,dum3,dum4,dum5
real, dimension(:,:,:,:), allocatable :: u,v,emc
real, dimension(:,:,:), allocatable :: mt
real, dimension(:,:,:), allocatable :: um,uma,vm,vma,emcm,emcma
real, dimension(:,:), allocatable :: mtma,mtm
real, dimension(:,:), allocatable :: ua,va
real, dimension(:,:,:), allocatable :: reg_mei
real, dimension(:), allocatable :: mei,mei_mon,tmask
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

tnt=0
do yr=yrs,yre
  mo_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  if (mod(yr,4)==0) mo_da(2)=mo_da(2)+1
  tnt=tnt+sum(mo_da,1)
enddo
write(*,*) tnt

! create mask
allocate(tmask(tnt))
tmask=-1
n=1
do yr=yrs,yre
  mo_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  if (mod(yr,4)==0) mo_da(2)=mo_da(2)+1
  nt=sum(mo_da,1)

  if (yr==yrs) then
    i=n-1+sum(mo_da(1:11),1)+1-10
    j=n-1+nt
    tmask(i:j)=1 
  elseif (yr==yre) then
    i=n-1+1
    j=n-1+sum(mo_da(1:3),1)+10
    tmask(i:j)=1
  else
    i=n-1+1
    j=n-1+sum(mo_da(1:3),1)+10
    tmask(i:j)=1
    i=n-1+sum(mo_da(1:11),1)+1-10
    j=n-1+nt
    tmask(i:j)=1 
  endif

  if (mod(yr,4)==0) then
    i=n-1+sum(mo_da(1:2),1)
    tmask(i)=-1
  endif

  n=n+nt
enddo

!write(*,*) sum(tmask,mask=tmask>0)

i=minloc(abs(lat-lats),1)
j=minloc(abs(lat-late),1)
sny=j-i+1

start=(/1,i,1/)

allocate(um(sny,nz,365),vm(sny,nz,365),emcm(sny,nz,365),mtm(sny,365))

fname=trim(path)//"/data/ERA5_spectrum_mean.dat"
open(10,file=trim(fname),access="direct",recl=sny*nz*3+sny)
do t=1,365
  read(10,rec=t) um(:,:,t),vm(:,:,t),emcm(:,:,t),mtm(:,t)
enddo
close(10)

allocate(uma(sny,nz,tnt),emcma(sny,nz,tnt),mtma(sny,tnt))
allocate(ua(sny,nz),va(sny,nz))

uma=0.
emcma=0.
mtma=0.

n=1
do yr=yrs,yre
  mo_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  if (mod(yr,4)==0) mo_da(2)=mo_da(2)+1
  nt=sum(mo_da,1)
  count=(/nx,sny,nt/)
  allocate(u(nx,sny,nz,nt),v(nx,sny,nz,nt),emc(nx,sny,nz,nt),mt(nx,sny,nt))

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

    write(fname,'(4A,I4,A)') trim(path),"/ERA5/v/v",trim(adjustl(tdum1)),"_",yr,".nc"
    !write(*,*) trim(fname),access(trim(fname),' ')

    ierr=nf90_open(trim(fname),nf90_nowrite,ncid1)
    if (ierr/=nf90_noerr) write(*,*) "open fail"
    ierr=nf90_inq_varid(ncid1,'v',varid1)
    if (ierr/=nf90_noerr) write(*,*) "inq var fail"
    ierr=nf90_get_var(ncid1,varid1,v(:,:,k,:),start=start,count=count)
    if (ierr/=nf90_noerr) write(*,*) "read fail"
    ierr=nf90_close(ncid1)
    if (ierr/=nf90_noerr) write(*,*) "close fail"
  enddo

  write(fname,'(2A,I4,A)') trim(path),"/../../ERA5/MT/MT_",yr,".nc"
 
  ierr=nf90_open(trim(fname),nf90_nowrite,ncid1)
  if (ierr/=nf90_noerr) write(*,*) "open fail"
  ierr=nf90_inq_varid(ncid1,'mt',varid1)
  if (ierr/=nf90_noerr) write(*,*) "inq var fail"
  ierr=nf90_get_var(ncid1,varid1,mt(:,:,:),start=start,count=count)
  if (ierr/=nf90_noerr) write(*,*) "read fail"
  ierr=nf90_close(ncid1)
  if (ierr/=nf90_noerr) write(*,*) "close fail" 
 

  do t=1,nt
    ua(:,:)=sum(u(:,:,:,t),1)/real(nx)
    va(:,:)=sum(v(:,:,:,t),1)/real(nx)
    do i=1,nx
      emc(i,:,:,t)=(u(i,:,:,t)-ua(:,:))*(v(i,:,:,t)-va(:,:))
    enddo
 
    o=n+t-1
    mtma(:,o)=sum(mt(:,:,t),1)/real(nx)
    uma(:,:,o)=sum(u(:,:,:,t),1)/real(nx)
    emcma(:,:,o)=sum(emc(:,:,:,t),1)/real(nx)
  enddo

  n=n+nt
  deallocate(u,v,emc,mt)
  write(*,*) yr
enddo

!goto 200

! linearly remove ENSO effect

allocate(mei_mon((yre-yrs+1)*12))
fname=trim(path)//"/MEI_v2.txt"
open(10,file=trim(fname))
read(10,*) 
i=1
do yr=yrs,yre
  read(10,'(I4,12F9.2)') j, mei_mon(i:i+12-1) 
  i=i+12
enddo
close(10)
!write(*,*) mei

! interplate MEI index to daily resolution
allocate(mei(tnt))
n=1
m=1
mei=-999.
do yr=yrs,yre
  mo_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  if (mod(yr,4)==0) mo_da(2)=mo_da(2)+1
  nt=sum(mo_da,1)

  do i=1,12
    o=n-1+sum(mo_da(1:i),1)-mo_da(i)+15
    j=m-1+i
    mei(o)=mei_mon(j)
    !if (yr==2021 .and. i==12) write(*,*) o
  enddo
  n=n+nt
  m=m+12
enddo

!write(*,*) "month",mei(1:60)
i=0
j=0
do t=1,tnt
  if (mei(t)>-100.) then
    if (i==0) then
      i=t
      mei(1:i-1)=mei(i)
      j=i
    else
      i=t
      do k=j+1,i-1
        mei(k)=mei(j)+(mei(i)-mei(j))/real(i-j+1)*real(k-j+1)
      enddo
      j=i
    endif
  endif
enddo
!write(*,*) "inter",mei(1:60)

allocate(reg_mei(sny,nz,3))
fname=trim(path)//"/data/ERA5_month_reg_mei.dat"
open(10,file=trim(fname),access="direct",recl=sny*nz*2+sny)
read(10,rec=1) reg_mei(:,:,1),reg_mei(:,:,2),reg_mei(:,1,3)
close(10)

do t=1,tnt
  uma(:,:,t)=uma(:,:,t)-mei(t)*reg_mei(:,:,1) 
  emcma(:,:,t)=emcma(:,:,t)-mei(t)*reg_mei(:,:,2) 
  mtma(:,t)=mtma(:,t)-mei(t)*reg_mei(:,1,3) 
enddo

!200 continue

fname=trim(path)//"/data/ERA5_daily_mean.dat"
open(10,file=trim(fname),access="direct",recl=sny*nz*2+sny)
n=1
do t=1,tnt
  if (tmask(t)>0) then
    write(10,rec=n) uma(:,:,t),emcma(:,:,t),mtma(:,t)
    n=n+1
  endif
enddo
close(10)
write(*,*) n-1

end program daily_zonal_mean
