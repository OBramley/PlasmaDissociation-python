program dynamics

implicit none
integer :: nst, inist, natom, ndim, n1, n2, i, j, k,n, nr, im, branch,ist,nb,natoms,skip
real (kind=8), allocatable, dimension (:) :: ElPhase
real (kind=8) CompTime, TimeStep, time, E0, E1, Mau,  val,prevtime,tocomplex
complex(8)    :: ii
real (kind=8), allocatable, dimension (:) :: R0,R1,P0,P1,V0,V1,F0,F1,Fb,M,Es0,Es1,Esb,tb
real (kind=8), allocatable, dimension (:,:) :: Fs0,Fs1,Fsb
real (kind=8), allocatable, dimension (:,:,:) :: Cs0,Cs1,Csb
complex (kind=8), allocatable, dimension (:) :: A0,A1,Ab
complex (kind=8), allocatable, dimension (:,:) :: HE_0,HE_1,HE_b
character (len=2), allocatable, dimension (:) :: atom
character (len=50) f_inp
ii=(0.,1.)
Mau=1822.887
call getarg(1,f_inp)

open(1,file=trim(f_inp)//'1.all')
read(1,*) natom, nst
read(1,*) branch
read(1,*) time, TimeStep
ndim=natom*3
prevtime=time

!Allocation
allocate(r0(ndim),p0(ndim),r1(ndim),p1(ndim), atom(natom), tb(natom))



tb=0.1e20
open(2,file=trim(f_inp)//'.diss')
do
  read(2,*,iostat=ist) val,nb
  if (ist==0) then
     do i=1,natom
       if (tb(i).lt.0.1e19.and.nb.ge.i) nb=nb+1
     end do
     tb(nb)=val
  else
     exit
  end if
end do
close(2)




! Read geometry
read(1,*)
do i=1,natom
  read(1,*) atom(i),r0(i*3-2:i*3)
end do
read(1,*)
do i=1,natom
  read(1,*) p0(i*3-2:i*3)
end do

!Allocation
allocate (A0(nst), A1(nst), Ab(nst))
allocate (Es0(nst),Fs0(ndim,nst))
allocate (Cs0(ndim,nst,nst))
allocate (HE_0(nst,nst))
allocate (F0(ndim),V0(ndim),V1(ndim))
allocate (M(ndim))

! Read amplitudes
read(1,*)
do i=1,nst
  read(1,*,iostat=ist) A0(i)
  if(ist.ne.0)then 
    backspace(1)
    read(1,*) tocomplex
    A0(i) = cmplx(tocomplex,0.0)
    end if 
end do

! Set masses 
do i=1,natom
  if (atom(i)=='C') then 
    M(i*3-2:i*3)=12*Mau
  else if (atom(i)=='Si') then 
     M(i*3-2:i*3)=28*Mau
   else if (atom(i)=='O') then 
     M(i*3-2:i*3)=16*Mau
  else if (atom(i)=='N') then 
    M(i*3-2:i*3)=14*Mau
  else if (atom(i)=='H') then 
    M(i*3-2:i*3)=Mau
  else if (atom(i)=='D') then 
    M(i*3-2:i*3)=2*Mau
  else if (atom(i)=='F') then
    M(i*3-2:i*3)=19*Mau
  else
    Print*,' Atom ',atom(i), ' is not suppoted'
    stop
  end if
end do

! Read potential energies
read(1,*)
do n=1,nst
  read(1,*) val, i  
  Es0(i)=val
end do

! Read forces
read(1,*)
nr=nst*ndim
do n=1,nr
    read(1,*) val, i,j
    Fs0(j,i)=val
end do

! Read couplings
read(1,*)
nr=ndim*nst*(nst-1)/2
do n=1,nr

      read(1,*) val,i,k,j
      Cs0(j,i,k)=val
      Cs0(j,k,i)=-Cs0(j,i,k)

end do



read(1,*)



open(2,file=trim(f_inp)//'.pop')
write(2,'(f8.2,9f14.8)')time,abs(A0)**2


open(3,file=trim(f_inp)//'.enr')
write(3,'(f8.2,12e20.8)')time,0.5*sum(P0**2/M),sum(abs(A0)**2*Es0),Es0

open(4,file=trim(f_inp)//'.xyz')
write(4,'(i4)') natom
write(4,'(f10.2)') Time
do i=1,natom
  write(4,'(a2,3f20.12)') atom(i),r0(i*3-2:i*3)/1.88973
end do


do 

read(1,*,iostat=ist) natoms, nst
if(ist.ne.0) exit
read(1,*) branch
read(1,*) time, TimeStep

! Read geometry
read(1,*)
do i=1,natom
  if (time.lt.tb(i)-0.1*TimeStep) then
    read(1,*) atom(i),r0(i*3-2:i*3)
  else
   r0(i*3-2:i*3)=r0(i*3-2:i*3)+TimeStep*p0(i*3-2:i*3)/M(i*3-2:i*3)
  end if    
end do

read(1,*)
do i=1,natom
  if (time.lt.tb(i)) read(1,*) p0(i*3-2:i*3)
end do

! Read amplitudes
read(1,*)
do i=1,nst
  read(1,*) A0(i)
end do

read(1,*)
do n=1,nst
  read(1,*) val, i
  Es0(i)=val
end do

! Read forces
read(1,*)

nr=3*natoms*nst
do n=1,nr
    read(1,*) val, i,j
    Fs0(j,i)=val
end do

! Read couplings
read(1,*)
nr=3*natoms*nst*(nst-1)/2

do n=1,nr

      read(1,*) val,i,k,j
      Cs0(j,i,k)=val
      Cs0(j,k,i)=-Cs0(j,i,k)

end do
read(1,*)
if(time.eq.prevtime)then
  continue
end if
prevtime =time 

write(2,'(f8.2,9f14.8)')time,abs(A0)**2
write(3,'(f8.2,12e20.8)')time,0.5*sum(P0**2/M),sum(abs(A0)**2*Es0),Es0
write(4,'(i4)') natom
write(4,'(f10.2)') Time 
do i=1,natom
  write(4,'(a2,3f20.12)') atom(i),r0(i*3-2:i*3)/1.88973
end do



end do

close(1)
close(2)
close(3)
end


