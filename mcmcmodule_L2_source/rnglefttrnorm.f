c23456

C ===============================
C GEORGIOS KARAGIANNIS @BRISTOL.AC.UK
C GEORGIOS.KARAGIANNIS@PNNL.GOV
C PNNL
C 2010-01-01
C 2013-08-20
C ===============================


      subroutine rngLeftTrNorm(rnd , mu , sig , lower )

        implicit none

        double precision, intent(in)    :: mu
        double precision, intent(in)    :: sig
        double precision, intent(in)    :: lower
        double precision, intent(out)   :: rnd

Cf2py intent(in)    mu
Cf2py intent(in)    sig
Cf2py intent(in)    lower
Cf2py intent(out)   rnd

        double precision :: z_lower
        double precision :: ell_opt
        double precision :: rho
        double precision :: u

        ! GENERATE Z ~ TN(0,1) ;  Z > (LOWER-MU)/SIG

        z_lower = (lower - mu)/sig
        ell_opt = 0.5*( z_lower + sqrt(z_lower**2 + 4.0) )

        do
            call rngLeftTrExp(rnd, ell_opt, z_lower)
            rho = exp( -0.5*( rnd - ell_opt )**2 )
            call rnguniform(u)
            if (u .lt. rho) then
                exit
            end if
        end do

        ! GENERATE X ~ TN(MU,SIG) ; X > LOWER

        rnd = mu + sig*rnd

      end


