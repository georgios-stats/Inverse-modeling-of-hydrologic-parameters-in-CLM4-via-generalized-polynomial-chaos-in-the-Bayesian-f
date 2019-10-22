C ===============================
C GEORGIOS KARAGIANNIS @BRISTOL.AC.UK
C GEORGIOS.KARAGIANNIS@PNNL.GOV
C PNNL
C 2010-01-01
C 2013-08-20
C ===============================

      subroutine rngLeftTrExp(rnd,  ell , lower )

        implicit none

        double precision, intent(in)    :: lower
        double precision, intent(in)    :: ell
        double precision, intent(out)   :: rnd

Cf2py intent(in)    lower
Cf2py intent(in)    ell
Cf2py intent(out)   rnd

        call rnguniform(rnd)

        rnd = lower -1.0/ell *log(1.0-rnd)

      end subroutine rngLeftTrExp
