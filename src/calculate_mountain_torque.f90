program cal_mt
use mpi
use netcdf
implicit none 

integer, parameter :: nx=480,ny=241,nz=6
real, parameter :: grav=9.80665, p0=100000., a=6378000.
real, parameter :: tri_pi=4.*atan(1.), d2r=tri_pi/180.
integer :: nproc,myid
real, dimension(ny) :: lat
real, dimension(nx) :: lon
integer :: i,j,k,m,n,o,t
integer :: nt,tn,ts,te
integer :: ierr,ncid1,varid1,dimid1,ncid2,varid2,dimid2
integer :: nDimensions,nVariables,nAttributes,unlimitedDimId
integer :: odim4d(4)
integer, dimension(:), allocatable :: dim4d
integer :: yr,da,access
integer :: yrs,yre
integer :: is,ie
integer, dimension(3) :: start,count
integer, dimension(12) :: mo_da
integer, dimension(:), allocatable :: mpi_s, mpi_n
real, dimension(:,:,:), allocatable :: sp,mt
real, dimension(:,:), allocatable :: time_bnds
real, dimension(nx,ny) :: topo,dtopo
real :: lats,late,dlon 
real :: dum1,dum2,dum3,dum4,dum5
character(200) :: path,fname,tdum1,tdum2

call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,nproc,ierr)
call mpi_comm_rank(mpi_comm_world,myid,ierr)

allocate(mpi_s(nproc),mpi_n(nproc))

path="/data/der0318/work/am_NH/"

yrs=1979
yre=2017

nt=yre-yrs+1

tn=nt/nproc
i=mod(nt,nproc)
ts=myid*tn+yrs
if (myid<i) then
  ts=ts+myid
  tn=tn+1
else
  ts=ts+i
endif
te=ts+tn-1
tn=te-ts+1
write(*,*) myid,ts,te

write(fname,'(2A)') trim(path),"/ITM/surface_geo.nc"

ierr=nf90_open(trim(fname),nf90_nowrite,ncid1)
if (ierr/=nf90_noerr) write(*,*) "open fail"
ierr=nf90_inq_varid(ncid1,'latitude',varid1)
if (ierr/=nf90_noerr) write(*,*) "inq var fail"
ierr=nf90_get_var(ncid1,varid1,lat)
if (ierr/=nf90_noerr) write(*,*) "read fail"
ierr=nf90_inq_varid(ncid1,'longitude',varid1)
if (ierr/=nf90_noerr) write(*,*) "inq var fail"
ierr=nf90_get_var(ncid1,varid1,lon)
if (ierr/=nf90_noerr) write(*,*) "read fail"
ierr=nf90_inq_varid(ncid1,'z',varid1)
if (ierr/=nf90_noerr) write(*,*) "inq var fail"
ierr=nf90_get_var(ncid1,varid1,topo)
if (ierr/=nf90_noerr) write(*,*) "read fail"
ierr=nf90_get_att(ncid1, varid1, "scale_factor", dum1)
if (ierr/=nf90_noerr) write(*,*) "scale fail"
ierr=nf90_get_att(ncid1, varid1, "add_offset", dum2)
if (ierr/=nf90_noerr) write(*,*) "offset fail"
ierr=nf90_close(ncid1)
if (ierr/=nf90_noerr) write(*,*) "close fail"
topo=topo*dum1+dum2
topo=topo/grav

dlon=2*(lon(2)-lon(1))*d2r
do j=1,ny
do i=1,nx
  is=i-1
  ie=i+1
  if (is<1) is=is+nx
  if (ie>nx) ie=ie-nx
  dtopo(i,j)=(topo(ie,j)-topo(is,j))/(dlon*cos(lat(j)*d2r))
enddo
enddo

do yr=ts,te
  mo_da=(/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
  if (mod(yr,4)==0) mo_da(2)=mo_da(2)+1
  nt=sum(mo_da,1)
  allocate(sp(nx,ny,nt),mt(nx,ny,nt))
  
  write(fname,'(2A,I4,A)') trim(path),"/ERA-I/SP/daily_interim_SP_",yr,".nc"
  
  ierr=nf90_open(trim(fname),nf90_nowrite,ncid1)
  if (ierr/=nf90_noerr) write(*,*) "open fail"
  ierr=nf90_inq_varid(ncid1,'sp',varid1)
  if (ierr/=nf90_noerr) write(*,*) "inq var fail"
  ierr=nf90_get_var(ncid1,varid1,sp)
  if (ierr/=nf90_noerr) write(*,*) "read fail"

  do t=1,nt
    mt(:,:,t)=grav/p0*sp(:,:,t)*dtopo(:,:)/a
  enddo

  write(fname,'(2A,I4,A)') trim(path),"/ITM/MT/MT_",yr,".nc"

  ierr=nf90_create(trim(fname),&
         cmode=nf90_64bit_offset,ncid=ncid2)
  if (ierr/=nf90_noerr) write(*,*) "create fail"
  ierr=nf90_inquire(ncid1,nDimensions,nVariables,nAttributes,unlimitedDimId)
  if (ierr/=nf90_noerr) write(*,*) "inquire fail"
   
  do i=1,nDimensions
    ierr=nf90_inquire_dimension(ncid1,i,tdum1,j)
    if (ierr/=nf90_noerr) write(*,*) "inq dim ",i," fail"
    ierr=nf90_def_dim(ncid2,trim(tdum1),j,dimid1)
    if (ierr/=nf90_noerr) write(*,*) "def dim ",i," fail"
  enddo

  do i=1,nVariables
    n=0
    ierr=nf90_inquire_variable(ncid1,i,name=tdum1,nDims=nDimensions,dimids=odim4d,nAtts=nAttributes)

    allocate(dim4d(nDimensions))
    do j=1,nDimensions
      dim4d(j)=odim4d(j) 
    enddo
   
    if (trim(tdum1)=="sp") then
      tdum1="mt"
      n=1
    endif
    ierr=nf90_def_var(ncid2,trim(tdum1),nf90_real,dim4d,varid1)
    if (ierr/=nf90_noerr) write(*,*) "def var ",i," fail"
    
    do j=1,nAttributes
      ierr=nf90_inq_attname(ncid1,i,j,name=tdum1)
      if (ierr/=nf90_noerr) write(*,*) "inq att ",j," fail"
      if (n==0) then
        ierr=nf90_copy_att(ncid1,i,trim(tdum1),ncid2,varid1)
        if (ierr/=nf90_noerr) write(*,*) "copy att ",j," fail"
      else
        tdum2="aaa"
        if (tdum1=="standard_name") tdum2="mountain_torque"
        if (tdum1=="long_name") tdum2="Mountain torque"
        if (tdum1=="units") tdum2="m s**-2"
        if (tdum1=="_FillValue" .or. &
            tdum1=="add_offset" .or. &
            tdum1=="scale_factor") tdum2="bbb"

        if (tdum2=="aaa") then
          ierr=nf90_copy_att(ncid1,i,trim(tdum1),ncid2,varid1)
          if (ierr/=nf90_noerr) write(*,*) "copy att ",j," fail"
        elseif (tdum2=="bbb") then
        else
          ierr=nf90_put_att(ncid2,varid1,trim(tdum1),trim(tdum2))
          if (ierr/=nf90_noerr) write(*,*) "put att ",j," fail"
        endif
      endif
    enddo

    deallocate(dim4d)
  enddo

  ierr=nf90_enddef(ncid2)
  !write(*,*) "enddef", nf90_strerror(ierr)
 
  ierr=nf90_inq_varid(ncid2,"longitude",varid1) 
  if (ierr/=nf90_noerr) write(*,*) "inq lon fail" 
  ierr=nf90_put_var(ncid2,varid1,lon)
  if (ierr/=nf90_noerr) write(*,*) "put lon fail", nf90_strerror(ierr)
  
  ierr=nf90_inq_varid(ncid2,"latitude",varid1) 
  if (ierr/=nf90_noerr) write(*,*) "inq lat fail" 
  ierr=nf90_put_var(ncid2,varid1,lat)
  if (ierr/=nf90_noerr) write(*,*) "put lat fail" 
  
  !allocate(time_bnds(2,nt))

  !ierr=nf90_inq_varid(ncid1,"time_bnds",varid1)
  !if (ierr/=nf90_noerr) write(*,*) "inq tbnd fail"
  !ierr=nf90_get_var(ncid1,varid1,time_bnds)
  !if (ierr/=nf90_noerr) write(*,*) "read tbnd fail", nf90_strerror(ierr)
  !ierr=nf90_inq_varid(ncid2,"time_bnds",varid1)
  !if (ierr/=nf90_noerr) write(*,*) "inq tbnd fail"
  !ierr=nf90_put_var(ncid2,varid1,time_bnds)
  !if (ierr/=nf90_noerr) write(*,*) "put tbnd fail", nf90_strerror(ierr)

  !deallocate(time_bnds)

  allocate(time_bnds(1,nt))

  ierr=nf90_inq_varid(ncid1,"time",varid1)
  if (ierr/=nf90_noerr) write(*,*) "inq time fail"
  ierr=nf90_get_var(ncid1,varid1,time_bnds(1,:))
  if (ierr/=nf90_noerr) write(*,*) "read time fail", nf90_strerror(ierr)
  ierr=nf90_inq_varid(ncid2,"time",varid1)
  if (ierr/=nf90_noerr) write(*,*) "inq time fail"
  ierr=nf90_put_var(ncid2,varid1,time_bnds(1,:))
  if (ierr/=nf90_noerr) write(*,*) "put time fail", nf90_strerror(ierr)

  deallocate(time_bnds)

  ierr=nf90_inq_varid(ncid2,"mt",varid1)
  if (ierr/=nf90_noerr) write(*,*) "inq mt fail"
  ierr=nf90_put_var(ncid2,varid1,mt)
  if (ierr/=nf90_noerr) write(*,*) "put mt fail"

  ierr=nf90_close(ncid1)
  if (ierr/=nf90_noerr) write(*,*) "close fail 1" 
  
  ierr=nf90_close(ncid2)
  if (ierr/=nf90_noerr) write(*,*) "close fail 2" 

  deallocate(sp,mt)
  write(*,*) myid,yr
enddo


call mpi_finalize(ierr)

end program cal_mt
