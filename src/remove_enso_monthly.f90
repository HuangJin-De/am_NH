program monthly_zonal_mean
use netcdf
implicit none

integer, parameter :: nx=480,ny=241,nz=27
real, parameter :: tri_pi=4.*atan(1.)
real, dimension(nz) :: lev
real, dimension(ny) :: lat
real, dimension(nx) :: lon
integer :: i,j,k,m,n,o,t,ii,jj
integer :: nt,sny,snz,tnt
integer :: ierr,ncid1,varid1
integer :: yr,da,access
integer :: yrs,yre
integer, dimension(3) :: start3,count3
integer, dimension(4) :: start4,count4
integer, dimension(12) :: mo_da
real :: lats,late
real :: dum1,dum2,dum3,dum4,dum5,dum6,dum7,dum8
real, dimension(:,:,:,:), allocatable :: u,v,emc
real, dimension(:,:,:), allocatable :: um,uma,umw,vm,vmw,emcm,emcma,emcmw
real, dimension(:,:), allocatable :: ua,va
real, dimension(:,:,:), allocatable :: reg_mei
real, dimension(:,:,:), allocatable :: mt
real, dimension(:,:), allocatable :: mta,mtm,mtmw
real, dimension(:), allocatable :: mei
character(200) :: path,fname,tdum1

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

lev=lev*100

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


tnt=(yre-yrs+1)*12

allocate(um(sny,snz,365),umw(sny,snz,12))
allocate(vm(sny,snz,365),vmw(sny,snz,12))
allocate(emcm(sny,snz,365),emcmw(sny,snz,12))
allocate(mtm(sny,365),mtmw(sny,12))

fname=trim(path)//"/data/ERA5_spectrum_mean.dat"
open(10,file=trim(fname),access="direct",recl=sny*snz*3+sny)
do t=1,365
  read(10,rec=t) um(:,:,t),vm(:,:,t),emcm(:,:,t),mtm(:,t)
enddo
close(10)

mo_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
do i=1,12
  j=sum(mo_da(1:i))-mo_da(i)+1
  k=sum(mo_da(1:i))
  umw(:,:,i)=sum(um(:,:,j:k),3)/real(k-j+1)
  vmw(:,:,i)=sum(vm(:,:,j:k),3)/real(k-j+1)
  emcmw(:,:,i)=sum(emcm(:,:,j:k),3)/real(k-j+1)
  mtmw(:,i)=sum(mtm(:,j:k),2)/real(k-j+1)
enddo

allocate(uma(sny,snz,tnt),emcma(sny,snz,tnt),mta(sny,tnt))
allocate(ua(sny,snz),va(sny,snz))

uma=0.; emcma=0.; mta=0.

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

  write(fname,'(2A,I4,A)') trim(path),"/ITM/MT/MT_",yr,".nc"

  ierr=nf90_open(trim(fname),nf90_nowrite,ncid1)
  if (ierr/=nf90_noerr) write(*,*) "open fail"
  ierr=nf90_inq_varid(ncid1,'mt',varid1)
  if (ierr/=nf90_noerr) write(*,*) "inq var fail"
  ierr=nf90_get_var(ncid1,varid1,mt,start=start3,count=count3)
  if (ierr/=nf90_noerr) write(*,*) "read fail"
  ierr=nf90_close(ncid1)
  if (ierr/=nf90_noerr) write(*,*) "close fail"


  do t=1,nt
    ua(:,:)=sum(u(:,:,:,t),1)/real(nx)
    va(:,:)=sum(v(:,:,:,t),1)/real(nx)
    do i=1,nx
      emc(i,:,:,t)=(u(i,:,:,t)-ua(:,:))*(v(i,:,:,t)-va(:,:))
    enddo
  enddo

  do i=1,12
    j=sum(mo_da(1:i))-mo_da(i)+1
    k=sum(mo_da(1:i))
    o=n+i-1
    mta(:,o)=sum(sum(mt(:,:,j:k),3),1)/real((k-j+1)*nx)-mtmw(:,i)
    uma(:,:,o)=sum(sum(u(:,:,:,j:k),4),1)/real((k-j+1)*nx)-umw(:,:,i)
    emcma(:,:,o)=sum(sum(emc(:,:,:,j:k),4),1)/real((k-j+1)*nx)-emcmw(:,:,i)
  enddo

  n=n+12
  deallocate(u,v,emc,mt)
  write(*,*) yr
enddo


! linearly remove ENSO effect

allocate(mei(tnt))
write(*,*) shape(mei)
fname=trim(path)//"/MEI_v2.txt"
open(10,file=trim(fname))
read(10,*)
i=1
do yr=yrs,yre
  read(10,'(I4,12F9.2)') j, mei(i:i+12-1) 
  i=i+12
enddo
close(10)
!write(*,*) mei

! linear regression
allocate(reg_mei(sny,snz,3))
do k=1,snz
do j=1,sny
  dum3=0.
  dum4=0.
  dum5=0.
  m=0
  do i=1,tnt
    if (mod(i,12)<=3) then
      dum3=dum3+uma(j,k,i)
      dum4=dum4+mei(i)
      dum5=dum5+emcma(j,k,i)
      m=m+1
    endif   
  enddo
  dum3=dum3/real(m)
  dum4=dum4/real(m)
  dum1=0.
  dum2=0.
  dum6=0.
  do i=1,tnt
    if (mod(i,12)<=3) then
      dum1=dum1+(uma(j,k,i)-dum3)*(mei(i)-dum4)
      dum2=dum2+(mei(i)-dum4)**2.
      dum6=dum6+(emcma(j,k,i)-dum5)*(mei(i)-dum4)
    endif
  enddo
  reg_mei(j,k,1)=dum1/dum2
  reg_mei(j,k,2)=dum6/dum2
 
  if (k==1) then
    dum3=0.
    m=0
    do i=1,tnt
      if (mod(i,12)<=3) then
        dum3=dum3+mta(j,i)
        m=m+1
      endif
    enddo
    dum3=dum3/real(m)
    dum1=0.
    dum2=0.
    do i=1,tnt
      if (mod(i,12)<=3) then
        dum1=dum1+(mta(j,i)-dum3)*(mei(i)-dum4)
        dum2=dum2+(mei(i)-dum4)**2.
      endif
    enddo
    reg_mei(j,k,3)=dum1/dum2
  endif
enddo
enddo

do t=1,tnt
  if (mod(t,12)<=3) then
    uma(:,:,t)=uma(:,:,t)-mei(t)*reg_mei(:,:,1)
    emcma(:,:,t)=emcma(:,:,t)-mei(t)*reg_mei(:,:,2)
    mta(:,t)=mta(:,t)-mei(t)*reg_mei(:,1,3)
  endif
enddo


!fname=trim(path)//"/data/ERA5_month_spectrum_mean.dat"
!open(10,file=trim(fname),access="direct",recl=sny*snz)
!do i=1,12
!  write(10,rec=i) um(:,:,i)
!enddo
!close(10)

fname=trim(path)//"/data/ERA5_month_reg_mei.dat"
open(10,file=trim(fname),access="direct",recl=sny*snz*2+sny)
write(10,rec=1) reg_mei(:,:,1),reg_mei(:,:,2),reg_mei(:,1,3)
close(10)

fname=trim(path)//"/data/ERA5_month_mean.dat"
open(10,file=trim(fname),access="direct",recl=sny*snz)
n=1
do t=1,tnt
  if (mod(t,12)<=3) then
    write(10,rec=n) uma(:,:,t)
    n=n+1
  endif
enddo
close(10)
write(*,*) n

deallocate(um,uma,reg_mei)

end program monthly_zonal_mean
