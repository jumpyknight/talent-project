      call pmat
      stop 
      end

c construct pairing Hamiltonian matrix
      subroutine pmat

c GLOSSARY:
c______________________________________________________________________
c
c b2bd    : basis of all 2-p states
c b4bd    : 4-p basis of paired states
c
c i(w/o #): counters
c
c j       : counter
c
c ham1    : diagonal part of hamiltonian matrix
c hf      : hamiltonian strength array
c hos     : hamiltonian occupied state array
c
c mps     : array containing any 4-p states
c
c nl      : number of 1-p spin-degenerate energy levels
c ns      : total number of 1-p states
c
c spbl    : array of 1-p energy levels
c spbs    : array of spins
c
c______________________________________________________________________

      implicit none
      integer i,j,nl,i1,i2,i3,i4,ns
      parameter(nl=4,ns=2*nl)
      integer :: spbl(ns),spbs(ns),mps(ns),b4bd(4,6),b2bd(2,28)
      integer :: hos(8),hf(8),ham1(6,6)

      i=0

 1000 format('1p states:')
 1001 format('4p basis:')
 1002 format('1p states in 4p basis:')
 1003 format('2p basis:')
 1004 format('1p states in 2p basis:')
 1005 format(502i2)
 1006 format('1p hamiltonian strength & occupied level:')
 1007 format('1p hamiltonian matrix')

      open(unit=1,file='pmat.out',status='unknown')

c all 1-p quantum numbers
      do j=1,nl
         i=i+1
         spbl(i)=j-1
         spbs(i)=1
         i=i+1
         spbl(i)=j-1
         spbs(i)=-1
      enddo

c write 1-p quantum numbers
      write(1,*)
      write(*,*)
      write(1,1000)
      write(*,1000)
      write(1,*)
      write(*,*)
      do i=1,8
         write(1,*)spbl(i),spbs(i)
         write(*,*)spbl(i),spbs(i)
      enddo
      write(1,*)

      mps(:)=0

c find 4-p basis
      write(1,1001)
      write(*,1001)
      write(1,*)
      write(*,*)
      j=0
      b4bd(:,:)=0
      do i1 = 1 , ns
       do i2 = i1+1 , ns
        do i3 = i2+1 , ns
         do i4 = i3+1 , ns
          mps(:) = 0
          mps(i1) = 1
          mps(i2) = 1
          mps(i3) = 1
          mps(i4) = 1
          if(spbs(i1).ne.-1*spbs(i2)) goto 10
          if(spbl(i1).ne.spbl(i2)) goto 10
          if(spbs(i3).ne.-1*spbs(i4)) goto 10
          if(spbl(i3).ne.spbl(i4)) goto 10
          write(1,*)(mps(i),i=1,8)
          write(*,*)(mps(i),i=1,8)
          j=j+1
          b4bd(1,j)=i1
          b4bd(2,j)=i2
          b4bd(3,j)=i3
          b4bd(4,j)=i4
   10 continue
      enddo
      enddo
      enddo
      enddo
      write(1,*)
      write(*,*)

      write(1,1002)
      write(*,1002)
      write(1,*)
      write(*,*)
      do i=1,4
         write(1,*)(b4bd(i,j),j=1,6)
         write(*,*)(b4bd(i,j),j=1,6)
      enddo
      write(1,*)
      write(*,*)

      write(1,1003)
      write(*,1003)
      write(1,*)
      write(*,*)

c TWO-BODY CONFIGURATIONS
      j=0
      b2bd(:,:)=0
      do i1=1,ns
         do i2=i1+1,ns
            mps(:)=0
            mps(i1)=1
            mps(i2)=1
            if(spbl(i1).eq.spbl(i2).and.spbs(i1).eq.spbs(i2)) goto 20
            write(1,*)(mps(i),i=1,8)
            write(*,*)(mps(i),i=1,8)
            j=j+1
            b2bd(1,j)=i1
            b2bd(2,j)=i2
   20       continue
         enddo
      enddo

      write(1,*)
      write(*,*)
      write(1,1004)
      write(*,1004)
      write(1,*)
      write(*,*)

      do i=1,2
         write(*,*)
         write(*,1005)(b2bd(i,j),j=1,28)
         write(1,*)
         write(1,1005)(b2bd(i,j),j=1,28)
      enddo

      write(*,*)

c construct 1p hamiltonian
      do i=1,ns
         hos(i)=i
      enddo

      i=1

      do j=1,4
         hf(i)=j-1
         hf(i+1)=j-1
         i=i+2
      enddo

      write(*,*)
      write(*,1006)
      write(*,*)
      write(1,*)
      write(1,1006)
      write(1,*)
      do i=1,8
         write(*,*)hf(i),hos(i)
         write(1,*)hf(i),hos(i)
      enddo

c calculate 1p matrix elements
      call obme(b4bd,hf,hos,ham1)
      write(*,*)
      write(*,1007)
      write(*,*)
      write(1,*)
      write(1,1007)
      write(1,*)
      do i=1,6
         write(*,*)(ham1(i,j),j=1,6)
         write(1,*)(ham1(i,j),j=1,6)
      enddo

      close(1)

      return
      end

c 1body matrix elements
      subroutine obme(b4bd,hf,hos,ham1)

c GLOSSARY:
c______________________________________________________________________
c
c a       : counter
c
c b       : counter
c b4bd    : 4p basis matrix
c 
c ham1    : diagonal part of the hamiltonian matrix
c hf      : hamiltonian strength array
c hos     : hamiltonian occupied state array
c
c i       : counter
c
c j       : counter
c
c k       : counter
c
c l       : counter
c
c m       : counter
c mel     : hamiltonian matrix element
c melt    : temporary hamiltonian matrix element
c______________________________________________________________________

      implicit none

      integer :: i,j,k,l,m,a,b
      integer :: b4bd(4,6),hf(8),hos(8),ham1(6,6)
      real    :: mel,melt


c i counters bras
c j counters ket
      do i=1,6					
         do j=1,6
            mel=0.d0
            do l=1,4
               do m=1,4
                  if(b4bd(l,i).eq.b4bd(m,j))melt=hf(b4bd(l,i))
                  if(b4bd(l,i).ne.b4bd(m,j))melt=0.
                  a=1
                  b=1
                  do k=1,3
                    if(a.eq.l)a=a+1
                    if(b.eq.m)b=b+1
                    if(b4bd(a,i).ne.b4bd(b,j))melt=0.
                    a=a+1
                    b=b+1
                  enddo
                  mel=mel+melt
               enddo
            enddo
         ham1(i,j)=mel
         enddo
      enddo

      return
      end
