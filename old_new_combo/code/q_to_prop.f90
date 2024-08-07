implicit none
character l2*63,l2t*63,l1t*38,l1*60, ic*1
integer i, j, k, natom, nst, ndim, iat, iend, ia
real(kind=8), allocatable ::  c(:,:,:), f(:,:), e(:)



open(1,file='t.0')
read(1,*) natom, nst

if (nst>9) then
  print*, 'ERROR: nimber of states >9'
  stop
end if

ndim=3*natom
allocate(c(ndim,nst,nst), f(ndim,nst), e(nst))

C=0.

l1t=' Total energy in the final basis set = '



l2t=' Gradient of SCF Energy'

!do i=1,nst
!write(ic,'(i1)') i

!open(11,file='f'//ic//'.out')

open(11,file='f.out')

f=0.
e=0.

do
  read(11,'(a60)')l1
  if(l1(1:38).eq.l1t) then
      read(l1(39:),*) e(1)
    exit
  end if
enddo

do
  read(11,'(a23)')l2
  if(l2.eq.l2t) then
    do j=1,(natom-1)/6+1
      read(11,*)
      iend=min(j*6,natom)
      read(11,*) k,f(18*(j-1)+1:3*iend:3,1)
      read(11,*) k,f(18*(j-1)+2:3*iend:3,1)
      read(11,*) k,f(18*(j-1)+3:3*iend:3,1)
    end do
    exit
  end if
enddo
close(11)



f=-f


!Write electronic structure data
do i=1,nst
  print'(e25.16,i8)', E(i), i
end do

print*
do i=1,nst
  do j=1,ndim
    print'(e25.16,2i8)', F(j,i), i,j
  end do
end do

print*
do i=1,nst-1
  do k= i+1,nst
    do j=1,ndim
      print'(e25.16,3i8)', C(j,i,k), i,k,j
    end do
  end do
end do

end

