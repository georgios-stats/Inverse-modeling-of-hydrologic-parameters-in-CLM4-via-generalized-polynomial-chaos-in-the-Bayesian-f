C MODIFIED FROM RANDLIB
C
C CORRECTED IN 2009/12/11
C
C c = a + int((b-a+1)*u)
C
C c ~DU(A,B), c \in \{a,...,b\}
C 

C ===============================
C GEORGIOS KARAGIANNIS @BRISTOL.AC.UK
C GEORGIOS.KARAGIANNIS@PNNL.GOV
C PNNL
C 2010-01-01
C 2013-08-20
C ===============================


      subroutine genprm(ivec,n)
C**********************************************************************
C
C    SUBROUTINE rnguniform( ivec, n )
C               GENerate random PeRMutation of ivec
C
C
C                              Arguments
C
C
C     ixivecOn output ivec is a random permutation of its
C                 value on input
C                         INTEGER ivec( n )
C
C     n <--> Length of ivec
C                         INTEGER n
C
C**********************************************************************
C     .. Scalar Arguments ..
      INTEGER n
C     ..
C     .. Array Arguments ..
      INTEGER ivec(n)
C     ..
C
Cf2py intent(in)  n 
Cf2py intent(in,out) ivec
C
C     .. Local Scalars ..
      INTEGER i,itmp,j
C     ..
      DOUBLE PRECISION u
C     ..
C     .. Executable Statements ..
      do i = 1,n-1
          call rnguniform(u)
          j = i + int((n-i+1)*u)
          itmp = ivec(j)
          ivec(j) = ivec(i)
          ivec(i) = itmp
      enddo

      end subroutine





