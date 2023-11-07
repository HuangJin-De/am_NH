program pca_analysis
use netcdf
use forpca, only: tpca
use kinds, only: rk
implicit none

real(rk), dimension(:,:), allocatable :: matrix,matrix_app,coeff,score
real(rk), dimension(:),   allocatable :: latent, explained
type(tpca) :: p
integer, parameter :: nx=576,ny=360,nz=6
real, parameter :: tri_pi=4.*atan(1.),grav=9.8
real, dimension(nz), parameter :: lev=(/1000.,925.,850.,700.,500.,250./)
real, dimension(ny) :: lat
real, dimension(:), allocatable :: slat
real, dimension(nx) :: lon
integer :: i,j,k,m,n,o,t
integer :: nt,sny,tnt
integer :: ierr,ncid1,varid1
integer :: yr,da,access
integer :: yrs,yre
integer, dimension(3) :: start,count
integer, dimension(12) :: mo_da
real :: lats,late
real :: dum1,dum2,dum3,dum4,dum5
real, dimension(:,:,:), allocatable :: ua, um
real, dimension(:,:), allocatable :: intua, intum, intp
real, dimension(:,:,:), allocatable :: u_reg
real, dimension(:,:), allocatable :: pc 
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
!write(*,*) tnt

i=minloc(abs(lat-lats),1)
j=minloc(abs(lat-late),1)
sny=j-i+1

allocate(slat(sny))
allocate(ua(sny,nz,tnt),um(sny,nz,4))

slat=lat(i:j)


fname=trim(path)//"/data/ERA5_month_mean.dat"
open(10,file=trim(fname),access="direct",recl=sny*nz)
do t=1,tnt
  read(10,rec=t) ua(:,:,t)
enddo
close(10)

fname=trim(path)//"/data/ERA5_month_spectrum_mean.dat"
open(10,file=trim(fname),access="direct",recl=sny*nz)
n=1
do t=1,12
  if (mod(t,12)<=3) then
    read(10,rec=t) um(:,:,n)
    n=n+1
  endif
enddo
close(10)

! integrated zonal wind
allocate(intua(sny,tnt),intp(sny,tnt),intum(sny,4))
intua=0.
intum=0.
intp=0.
do t=1,tnt
  i=mod(t,4)
  if (i==0) i=4
  do j=1,sny
    do k=1,nz-1
      intua(j,t)=intua(j,t)+0.5*(ua(j,k  ,t)+um(j,k  ,i) &
                                +ua(j,k+1,t)+um(j,k+1,i))*(lev(k)-lev(k+1))/grav
      intp(j,t)=intp(j,t)+(lev(k)-lev(k+1))/grav
      if (t<=4) &
      intum(j,t)=intum(j,t)+0.5*(um(j,k  ,i) &
                                +um(j,k+1,i))*(lev(k)-lev(k+1))/grav
    enddo
    intua(j,t)=intua(j,t)-intum(j,i)
    intua(j,t)=intua(j,t)*cos(slat(j)*tri_pi/180.)/intp(j,t)
  enddo
enddo

fname=trim(path)//"/data/ERA5_intua.dat"
open(10,file=trim(fname),access="direct",recl=sny)
do t=1,tnt
  write(10,rec=t) intua(:,t)
enddo
close(10)

! perform PCA
allocate(matrix(sny,tnt))
matrix=reshape(intua,shape(matrix))
call p%pca(matrix,tnt,'eig',coeff,score,latent,explained,matrix_app)

write(*,*) "expla: ",sum(explained,1),explained(1:10)

! normalized PCs
do t=1,tnt
  dum1=sum(coeff(:,t),1)/real(tnt)
  dum2=0.
  do i=1,tnt
    dum2=dum2+(coeff(i,t)-dum1)**2.
  enddo
  dum2=sqrt(dum2/real(tnt))
  coeff(:,t)=(coeff(:,t)-dum1)/dum2
  score(:,t)=score(:,t)*dum2-dum1/dum2
enddo

! regress zonal wind to PC1
n=10
allocate(pc(tnt,n))
allocate(u_reg(sny,nz,n))
pc=coeff(:,1:n)
do m=1,n
do k=1,nz
do j=1,sny
  dum3=0.
  dum4=0.
  do i=1,tnt
    dum3=dum3+ua(j,k,i)
    dum4=dum4+pc(i,m)
  enddo
  dum3=dum3/real(tnt)
  dum4=dum4/real(tnt)
  dum1=0.
  dum2=0.
  do i=1,tnt
    dum1=dum1+(ua(j,k,i)-dum3)*(pc(i,m)-dum4)
    dum2=dum2+(pc(i,m)-dum4)**2.
  enddo
  u_reg(j,k,m)=dum1/dum2
enddo
enddo
enddo


fname=trim(path)//"/data/ERA5_month_reg_pc.dat"
open(10,file=trim(fname),access="direct",recl=sny*nz)
do m=1,n
  write(10,rec=m) u_reg(:,:,m)
enddo
close(10)

fname=trim(path)//"/data/ERA5_pc.dat"
open(10,file=trim(fname),access="direct",recl=sny)
do t=1,10
  write(10,rec=t) real(score(:,t))
enddo
close(10)

fname=trim(path)//"/data/ERA5_pc_series.dat"
open(10,file=trim(fname),access="direct",recl=tnt)
do t=1,tnt
  write(10,rec=t) real(coeff(:,t))
enddo
close(10)

call p%dlloc()
deallocate(ua,um,intua,intum)


end program pca_analysis
