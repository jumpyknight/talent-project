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

      i=0

 1000 format('1p states:')
 1001 format('4p basis:')
 1002 format('1p states in 4p basis:')
 1003 format('2p basis:')
 1004 format('1p states in 2p basis:')
 1005 format(502i2)

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
      enddo

      close(1)

      return
      end
      
      

      
      
      
      


