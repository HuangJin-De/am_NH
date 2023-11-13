program retopy
implicit none

integer :: i,j,k,m,n,o,t
integer :: nt,tnt
integer :: yr,da,access
integer :: yrs,yre,ts,te
real :: dum1,dum2,dum3,dum4,dum5,dum6
integer, dimension(12) :: mo_da
real, dimension(:,:), allocatable :: zk,mk
real, dimension(:,:), allocatable :: zp,mp
real, dimension(:,:), allocatable :: outvar
character(200) :: path,fname,tdum1

yrs=1979
yre=2021

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
allocate(zk(nt,2),mk(nt,2))
allocate(outvar(42,tnt-42*10))


n=1
do t=1,42
  ts=(t-1)*nt+1
  te=t*nt
    
  zk=zp(ts:te,1:2)
  mk=mp(ts:te,1:2)
   
  do i=11,121
    outvar(1,n)=mk(i,1)
    outvar(2,n)=mk(i,2)
  
    outvar(3:12,n)=zk(i-10:i-1,1)
    outvar(13:22,n)=zk(i-10:i-1,2)
    outvar(23:32,n)=mk(i-10:i-1,1)
    outvar(33:42,n)=mk(i-10:i-1,2)
   
    n=n+1
  enddo

enddo


fname="./train_data/zm_10days_data.dat"
open(10,file=trim(fname),access="direct",recl=42)
do i=1,tnt-10*42
  write(10,rec=i) outvar(:,i)
enddo
close(10)


! original data
fname="./data/ERA5_z_m_ori_data.dat"
open(10,file=trim(fname),access="direct",recl=20)
do t=1,tnt
  read(10,rec=t) zp(t,:),mp(t,:)
enddo
close(10)

nt=121

n=1
do t=1,42
  ts=(t-1)*nt+1
  te=t*nt

  zk=zp(ts:te,1:2)
  mk=mp(ts:te,1:2)

  do i=11,121
    outvar(1,n)=mk(i,1)
    outvar(2,n)=mk(i,2)

    outvar(3:12,n)=zk(i-10:i-1,1)
    outvar(13:22,n)=zk(i-10:i-1,2)
    outvar(23:32,n)=mk(i-10:i-1,1)
    outvar(33:42,n)=mk(i-10:i-1,2)

    n=n+1
  enddo

enddo


fname="./train_data/zm_10days_ori_data.dat"
open(10,file=trim(fname),access="direct",recl=42)
do i=1,tnt-10*42
  write(10,rec=i) outvar(:,i)
enddo
close(10)

end program retopy 
