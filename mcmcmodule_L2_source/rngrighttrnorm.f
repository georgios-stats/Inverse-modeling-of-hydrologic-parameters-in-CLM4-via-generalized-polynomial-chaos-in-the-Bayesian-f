C ===============================
C GEORGIOS KARAGIANNIS @BRISTOL.AC.UK
C GEORGIOS.KARAGIANNIS@PNNL.GOV
C PNNL
C 2010-01-01
C 2013-08-20
C ===============================

c23456
      subroutine rngRightTrNorm(rnd , mu , sig , upper )

        implicit none

        double precision, intent(in)    :: mu
        double precision, intent(in)    :: sig
        double precision, intent(in)    :: upper
        double precision, intent(out)   :: rnd

Cf2py intent(in)    mu
Cf2py intent(in)    sig
Cf2py intent(in)    lower
Cf2py intent(out)   rnd

        call rngLeftTrNorm( rnd, -mu , sig , -upper )

        rnd = -rnd

      end
