program projection_onto_eof
use netcdf
use, intrinsic :: iso_c_binding
implicit none

include 'fftw3.f03'

real, dimension(:,:), allocatable :: score
integer, parameter :: nx=576,ny=360,nz=6
real, parameter :: tri_pi=4.*atan(1.),grav=9.8,a=6378000.,d2r=tri_pi/180.
real, dimension(nz), parameter :: lev=(/1000.,925.,850.,700.,500.,250./)
real, dimension(ny) :: lat
real, dimension(:), allocatable :: slat
real, dimension(nx) :: lon
integer :: i,j,k,m,n,o,t
integer :: is,ie
integer :: nt,sny,tnt
integer :: ierr,ncid1,varid1
integer :: yr,da,access
integer :: yrs,yre
integer, dimension(3) :: start,count
integer, dimension(12) :: mo_da
real :: lats,late
real :: dum1,dum2,dum3,dum4,dum5
real, dimension(:), allocatable :: tmask,tjul,tyear,tmp1d
real, dimension(:), allocatable :: zk,mk,hann
real, dimension(:), allocatable :: zm,mm,zv,mv
real, dimension(:,:), allocatable :: zp,mp
real, dimension(:,:,:,:), allocatable :: csp
real, dimension(:,:,:), allocatable :: ua,emca,um,vm,emcm
real, dimension(:,:), allocatable :: mta,mtm
real, dimension(:,:), allocatable :: intua,intum,intemca,intemcm,intp
real, dimension(:,:,:), allocatable :: reg_mei
real, dimension(:,:), allocatable :: pc
character(200) :: path,fname,tdum1

! fft 
real(c_double), dimension(:), allocatable :: fin
complex(c_double_complex), dimension(:), allocatable :: fcoe, zn, mn
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

! create mask
allocate(tmask(tnt),tjul(tnt),tyear(tnt))
tmask=-1
tjul=-1
tyear=-1
n=1
do yr=yrs,yre
  mo_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  !if (mod(yr,4)==0) mo_da(2)=mo_da(2)+1
  nt=sum(mo_da,1)

  if (yr==yrs) then
    i=n
    j=n+sum(mo_da(12:12),1)
    tmask(i:j-1)=1
    tyear(i:j-1)=yr
    do k=sum(mo_da(1:11),1)+1,nt
      tjul(i+k-sum(mo_da(1:11),1)-1)=k
    enddo
    n=j
  elseif (yr==yre) then
    i=n
    j=n+sum(mo_da(1:3),1)
    tmask(i:j-1)=1
    tyear(i:j-1)=yr-1
    do k=1,sum(mo_da(1:3),1)
      tjul(i+k-1)=k
    enddo
    n=j
  else
    i=n
    j=n+sum(mo_da(1:3),1)
    tmask(i:j-1)=1
    tyear(i:j-1)=yr-1
    do k=1,sum(mo_da(1:3),1)
      tjul(i+k-1)=k
    enddo
    n=j

    i=n
    j=n+sum(mo_da(12:12),1)
    tmask(i:j-1)=1
    tyear(i:j-1)=yr
    do k=sum(mo_da(1:11),1)+1,nt
      tjul(i+k-sum(mo_da(1:11),1)-1)=k
    enddo
    n=j
  endif
enddo


allocate(um(sny,nz,365),vm(sny,nz,365),emcm(sny,nz,365),mtm(sny,365))

fname=trim(path)//"/data/ERA5_spectrum_mean.dat"
open(10,file=trim(fname),access="direct",recl=sny*nz*3+sny)
do t=1,365
  read(10,rec=t) um(:,:,t),vm(:,:,t),emcm(:,:,t),mtm(:,t)
enddo
close(10)

allocate(ua(sny,nz,tnt),emca(sny,nz,tnt),mta(sny,tnt))

fname=trim(path)//"/data/ERA5_daily_mean.dat"
open(10,file=trim(fname),access="direct",recl=sny*nz*2+sny)
do t=1,tnt
  read(10,rec=t) ua(:,:,t),emca(:,:,t),mta(:,t)
enddo
close(10)


allocate(score(sny,10))

fname=trim(path)//"/data/ERA5_pc.dat"
open(10,file=trim(fname),access="direct",recl=sny)
do t=1,10
  read(10,rec=t) score(:,t)
enddo
close(10)

allocate(intua(sny,tnt),intum(sny,tnt) &
        ,intemca(sny,tnt),intemcm(sny,tnt) &
        ,intp(sny,tnt))

intua=0.
intum=0.
intemca=0.
intemcm=0.

allocate(tmp1d(sny))
allocate(zp(tnt,10),mp(tnt,10))

m=1
do t=1,tnt
  i=tjul(t)
  if (i==-1) goto 989
  do j=1,sny
    do k=1,nz-1
      intua(j,t)=intua(j,t)+0.5*(ua(j,k,t)+ua(j,k+1,t))*(lev(k)-lev(k+1))/grav
      intemca(j,t)=intemca(j,t)+0.5*(emca(j,k,t)+emca(j,k+1,t))*(lev(k)-lev(k+1))/grav

      intp(j,t)=intp(j,t)+(lev(k)-lev(k+1))/grav

      intum(j,t)=intum(j,t)+0.5*(um(j,k,i)+um(j,k+1,i))*(lev(k)-lev(k+1))/grav
      intemcm(j,t)=intemcm(j,t)+0.5*(emcm(j,k,i)+emcm(j,k+1,i))*(lev(k)-lev(k+1))/grav
    enddo
  enddo
  
  tmp1d=0.
  tmp1d(1:sny-1)=0.5*(intemca(1:sny-1,t)*cos(slat(1:sny-1)*d2r)**2. &
                     +intemca(2:sny  ,t)*cos(slat(2:sny  )*d2r)**2.)
  intemca(2:sny-1,t)=(tmp1d(2:sny-1)-tmp1d(1:sny-2))/(slat(2:sny-1)-slat(1:sny-2))/a &
                      /cos(slat(2:sny-1)*d2r)**2
  intemca(1,t)=0.
  intemca(sny,t)=0.

  !write(*,*) "here:",yr,intemca(2:sny-1,t)

  tmp1d=0.
  tmp1d(1:sny-1)=0.5*(intemcm(1:sny-1,t)*cos(slat(1:sny-1)*d2r)**2. &
                     +intemcm(2:sny  ,t)*cos(slat(2:sny  )*d2r)**2.)
  intemcm(2:sny-1,t)=(tmp1d(2:sny-1)-tmp1d(1:sny-2))/(slat(2:sny-1)-slat(1:sny-2))/a &
                      /cos(slat(2:sny-1)*d2r)**2
  intemcm(1,t)=0.
  intemcm(sny,t)=0.

  intua(:,t)=intua(:,t)/intp(:,t)
  intum(:,t)=intum(:,t)/intp(:,t)
  intemca(:,t)=(intemca(:,t)+mta(:,t))/intp(:,t)
  intemcm(:,t)=(intemcm(:,t)+mtm(:,i))/intp(:,t)
 
  intemca(:,t)=-intemca(:,t)*86400.
  intemcm(:,t)=-intemcm(:,t)*86400.
 
  intua(:,t)=(intua(:,t)-intum(:,t))
  intemca(:,t)=(intemca(:,t)-intemcm(:,t))

  if (mod(t,121)==0) then
    open(10,file="./data/int_data.dat",access="direct",recl=sny*3*121)
    write(10,rec=m) intua(:,t-120:t)+intum(:,t-120:t) &
                   ,intemca(:,t-120:t)+intemcm(:,t-120:t) &
                   ,intp(:,t-120:t)
                 
    m=m+1
    close(10)
  endif
  
  !if (tyear(t)==2000 .and. (i==80 .or. i==81)) then 
  !  write(*,*) "ua:",intua(40,t),intum(40,t)
  !  write(*,*) "emca:",intemca(40,t),intemcm(40,t)
  !  write(*,*) "mta:",mta(40,t),mtm(40,i)
  !endif
 
  do n=1,10
    zp(t,n)=sum(intua(:,t)*cos(slat*d2r)*score(:,n))/sqrt(sum(score(:,n)*cos(slat*d2r)*score(:,n)))
    mp(t,n)=sum(intemca(:,t)*cos(slat*d2r)*score(:,n))/sqrt(sum(score(:,n)*cos(slat*d2r)*score(:,n)))
  enddo
  
  989 continue
enddo

!write(*,*) sqrt(sum(score(:,1)*cos(slat*d2r)*score(:,1)))

!allocate(zm(10),mm(10),zv(10),mv(10))
!zm=sum(zp,1)/real(tnt)
!mm=sum(mp,1)/real(tnt)
!
!do j=1,10
!  zv(j)=sqrt(sum((zp(:,j)-zm(j))**2.,1)/real(tnt))
!  zp(:,j)=(zp(:,j)-zm(j))/zv(j)
!
!  !if (j==1) write(*,*) zm(j),zv(j)
!
!  !mv(j)=sqrt(sum((mp(:,j)-zm(j))**2.,1)/real(tnt))
!  mp(:,j)=(mp(:,j)-zm(j))/zv(j)
!  !if (j==1) write(*,*) mm(j),mv(j)
!enddo

fname="./data/ERA5_z_m_ori_data.dat"
open(10,file=trim(fname),access="direct",recl=20)
do t=1,tnt
  write(10,rec=t) zp(t,:),mp(t,:)
enddo
close(10)

nt=121
deallocate(tmp1d)
allocate(tmp1d(nt))

do t=1,42
  is=(t-1)*nt+1
  ie=t*nt
  do j=1,10
    tmp1d(1:121)=zp(is:ie,j)
    dum1=sum(tmp1d(1:121),1)/real(nt)
    dum2=sqrt(sum((tmp1d(1:121)-dum1)**2.,1)/real(nt))
    zp(is:ie,j)=(tmp1d(1:121)-dum1)/dum2

    tmp1d(1:121)=mp(is:ie,j)
    dum1=sum(tmp1d(1:121),1)/real(nt)
    dum2=sqrt(sum((tmp1d(1:121)-dum1)**2.,1)/real(nt))
    mp(is:ie,j)=(tmp1d(1:121)-dum1)/dum2
    !if (j==10) write(*,*) dum1,dum2
  enddo
enddo 

fname="./data/ERA5_z_m_data.dat"
open(10,file=trim(fname),access="direct",recl=20)
do t=1,tnt
  write(10,rec=t) zp(t,:),mp(t,:)
enddo
close(10)


allocate(csp(61,10,6,42))
m=1
do yr=yrs,yre
  nt=sum(tyear,mask=((tyear==yr) .and. (tmask>0)))/yr
  !write(*,*) nt,nt/2+1
  if (nt<=0) goto 987
  allocate(zk(nt),mk(nt),hann(nt))
  allocate(fin(nt),fcoe(nt/2+1),zn(nt/2+1),mn(nt/2+1))

  dum1=real(nt+1)/2.
  do i=1,nt
    hann(i)=cos(tri_pi*(real(i)-dum1)/real(nt))**2.
    !if (yr==yrs) then
    !  write(*,*) i,hann(i)
    !endif
  enddo

  do i=1,10

    n=1
    do t=1,tnt
      if ((tyear(t)==yr) .and. (tmask(t)>0)) then
        zk(n)=zp(t,i)
        mk(n)=mp(t,i)
        n=n+1
        if (n>nt) goto 988 
      endif
    enddo
    988 continue

    !write(*,*) "zk:",zk
    !write(*,*) "mk:",mk

    plan_forward=fftw_plan_dft_r2c_1d(nt,fin,fcoe,FFTW_ESTIMATE)

    fin=zk*hann
    call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
    zn=fcoe
    
    fin=mk*hann
    call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
    mn=fcoe

    call fftw_free(plan_forward)

    if (i==1 .and. yr==yrs) then
      do j=1,1!nt/2+1
        write(*,*) j,real(j-1)/real(nt),sum(conjg(zn)*mn,1), &
                   sum(conjg(zn)*mn,1),sum(conjg(mn)*zn,1)!conjg(zn(j))*mn(j)/(conjg(zn(j))*zn(j))
      enddo
    endif

 
    csp(:,i,1,yr-yrs+1)=real(conjg(zn)*mn/(conjg(zn)*zn))
    csp(:,i,2,yr-yrs+1)=imag(conjg(zn)*mn/(conjg(zn)*zn))
    !csp(:,i,3,yr-yrs+1)=abs(zn)**2./sum(abs(zn)**2.,1)
    csp(:,i,3,yr-yrs+1)=abs(conjg(zn)*mn)**2/sum((conjg(zn)*zn)*(conjg(mn)*mn),1)
    !csp(:,i,4,yr-yrs+1)=abs(mn)**2./sum(abs(mn)**2.,1)
    csp(:,i,4,yr-yrs+1)=0.!abs(mn)**2./sum(abs(mn)**2.,1)

    plan_forward=fftw_plan_dft_r2c_1d(nt,fin,fcoe,FFTW_ESTIMATE)
    fin=zk
    call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
    zn=fcoe   
    call fftw_free(plan_forward)

    plan_forward=fftw_plan_dft_c2r_1d(nt,fcoe,fin,FFTW_ESTIMATE)
    fcoe=conjg(zn)*zn
    call fftw_execute_dft_c2r(plan_forward,fcoe,fin)
    call fftw_free(plan_forward)

    csp(:,i,5,yr-yrs+1)=fin(1:nt/2+1)/real(nt)**2

    plan_forward=fftw_plan_dft_r2c_1d(nt,fin,fcoe,FFTW_ESTIMATE)
    fin=mk
    call fftw_execute_dft_r2c(plan_forward,fin,fcoe)
    mn=fcoe   
    call fftw_free(plan_forward)

    plan_forward=fftw_plan_dft_c2r_1d(nt,fcoe,fin,FFTW_ESTIMATE)
    fcoe=conjg(mn)*mn
    call fftw_execute_dft_c2r(plan_forward,fcoe,fin)
    call fftw_free(plan_forward)

    csp(:,i,6,yr-yrs+1)=fin(1:nt/2+1)/real(nt)**2

  enddo

  open(10,file="./data/ERA5_crossspectrum_data.dat",access="direct",recl=6*(nt/2+1)*10)
  write(10,rec=m) csp(:,:,:,yr-yrs+1)
  m=m+1
  close(10)
 
  deallocate(zk,mk,hann,fin,fcoe,zn,mn)
enddo
987 continue 


end program projection_onto_eof
