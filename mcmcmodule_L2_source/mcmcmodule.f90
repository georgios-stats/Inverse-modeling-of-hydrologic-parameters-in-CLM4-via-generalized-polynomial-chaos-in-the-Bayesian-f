! ===============================
! GEORGIOS KARAGIANNIS (@PNNL.GOV)
! GEORGIOS.KARAGIANNIS@PNNL.GOV
! PNNL
! 2013-08-20
! ===============================


! DEPENDS ON:
!               rnguniform.f
!               rngnormal.f


! mcmcsampler =================================================================

subroutine mcmcsampler(Pr_ga, ga_sample, ga, &
                        rho_sample, rho, a_rho, b_rho, &
                        beta_sample, beta, &
                        sig2_sample, sig2, a_sig2, b_sig2, &
                        lambda_sample, lambda, a_lambda, b_lambda, &
                        dmax, y_en, YtY, XtY, XtX, &
                        Nsweep, Nburnin, &
                        fixed, &
                        sig2_thr &
                        )

    implicit none

    integer, intent(in)             :: dmax

    integer, intent(in)             :: Nsweep
    integer, intent(in)             :: Nburnin

    integer, intent(out)            :: ga_sample(dmax,Nsweep)
    integer, intent(inout)          :: ga(dmax)

    double precision, intent(out)   :: rho_sample(Nsweep)
    double precision, intent(inout) :: rho
    double precision, intent(in)    :: a_rho
    double precision, intent(in)    :: b_rho

    double precision, intent(out)   :: beta_sample(dmax,Nsweep)
    double precision, intent(inout) :: beta(dmax)

    double precision, intent(out)   :: sig2_sample(Nsweep)
    double precision, intent(inout) :: sig2
    double precision, intent(in)    :: a_sig2
    double precision, intent(in)    :: b_sig2

    double precision, intent(out)   :: lambda_sample(Nsweep)
    double precision, intent(inout) :: lambda
    double precision, intent(in)    :: a_lambda
    double precision, intent(in)    :: b_lambda

    integer, intent(in)             :: y_en
    double precision, intent(in)    :: YtY
    double precision, intent(in)    :: XtY(dmax)
    double precision, intent(in)    :: XtX(dmax,dmax)

    integer, intent(in)             :: fixed

    double precision, intent(in)    :: sig2_thr

    double precision, intent(out)   :: Pr_ga(dmax)

!f2py intent(in)        :: dmax
!f2py intent(in)        :: Nsweep
!f2py intent(in)        :: Nburnin
!f2py intent(out)       :: ga_sample
!f2py intent(in,out)    :: ga
!f2py intent(out)       :: rho_sample
!f2py intent(in,out)    :: rho
!f2py intent(in)        :: a_rho
!f2py intent(in)        :: b_rho
!f2py intent(out)       :: beta_sample
!f2py intent(in,out)    :: beta
!f2py intent(out)       :: sig2_sample
!f2py intent(in,out)    :: sig2
!f2py intent(in)        :: a_sig2
!f2py intent(in)        :: b_sig2
!f2py intent(out)       :: lambda_sample
!f2py intent(in,out)    :: lambda
!f2py intent(in)        :: a_lambda
!f2py intent(in)        :: b_lambda
!f2py intent(in)        :: y_en
!f2py intent(in)        :: YtY
!f2py intent(in)        :: XtY
!f2py intent(in)        :: XtX
!f2py intent(in)        :: fixed
!f2py intent(in)        :: sig2_thr
!f2py intent(out)       :: Pr_ga

    integer                         :: it
    integer                         :: j

    ! Burn in iterations

    do it = 1, Nburnin

        call mcmcsweep(Pr_ga(1:dmax), ga(1:dmax), &
                            rho, a_rho, b_rho, &
                            beta(1:dmax), &
                            sig2, a_sig2, b_sig2, &
                            lambda, a_lambda, b_lambda, &
                            dmax, y_en, YtY, XtY(1:dmax), XtX(1:dmax,1:dmax), &
                            fixed, &
                            sig2_thr &
                           )

    end do

    ! Sample iterations

    do it = 1, Nsweep

        call mcmcsweep(Pr_ga(1:dmax), ga(1:dmax), &
                            rho, a_rho, b_rho, &
                            beta(1:dmax), &
                            sig2, a_sig2, b_sig2, &
                            lambda, a_lambda, b_lambda, &
                            dmax, y_en, YtY, XtY(1:dmax), XtX(1:dmax,1:dmax), &
                            fixed, &
                            sig2_thr &
                           )

        do j = 1,dmax
            ga_sample(j,it) = ga(j)
            beta_sample(j,it) = beta(j)
        end do
        rho_sample(it) = rho
        sig2_sample(it) = sig2
        lambda_sample(it) = lambda

    end do

end subroutine mcmcsampler

! mcmcsweep ===================================================================

subroutine mcmcsweep(Pr_ga, ga, &
                        rho, a_rho, b_rho, &
                        beta, &
                        sig2, a_sig2, b_sig2, &
                        lambda, a_lambda, b_lambda, &
                        dmax, y_en, YtY, XtY, XtX, &
                        fixed, &
                        sig2_thr &
                       )

    implicit none

    integer, intent(in)             :: dmax

    integer, intent(inout)          :: ga(dmax)

    double precision, intent(inout) :: rho
    double precision, intent(in)    :: a_rho
    double precision, intent(in)    :: b_rho

    double precision, intent(inout) :: beta(dmax)

    double precision, intent(inout) :: sig2
    double precision, intent(in)    :: a_sig2
    double precision, intent(in)    :: b_sig2

    double precision, intent(inout) :: lambda
    double precision, intent(in)    :: a_lambda
    double precision, intent(in)    :: b_lambda

    integer, intent(in)             :: y_en
    double precision, intent(in)    :: YtY
    double precision, intent(in)    :: XtY(dmax)
    double precision, intent(in)    :: XtX(dmax,dmax)
    integer, intent(in)             :: fixed

    double precision, intent(in)    :: sig2_thr
    double precision, intent(out)   :: Pr_ga(dmax)

!f2py intent(in)        :: dmax
!f2py intent(in,out)    :: ga
!f2py intent(in,out)    :: rho
!f2py intent(in)        :: a_rho
!f2py intent(in)        :: b_rho
!f2py intent(in,out)    :: beta
!f2py intent(in,out)    :: sig2
!f2py intent(in)        :: a_sig2
!f2py intent(in)        :: b_sig2
!f2py intent(in,out)    :: lambda
!f2py intent(in)        :: a_lambda
!f2py intent(in)        :: b_lambda
!f2py intent(in)        :: y_en
!f2py intent(in)        :: YtY
!f2py intent(in)        :: XtY
!f2py intent(in)        :: XtX
!f2py intent(in)        :: sig2_thr
!f2py intent(out)       :: Pr_ga

    integer                         :: j
    integer                         :: i
    integer                         :: rep
    integer                         :: dga
    double precision                :: vrb
    double precision                :: vrbA
    double precision                :: vrbB
    double precision                :: sum_beta_sq
    double precision                :: BtXtXB
    double precision                :: BtXtY

    ! update ga, beta

do rep = 1, 10

    call update_ga_beta(Pr_ga(1:dmax), ga(1:dmax), &
                        rho, &
                        beta(1:dmax), &
                        sig2, &
                        lambda, &
                        dmax, XtY(1:dmax), XtX(1:dmax,1:dmax), &
                        fixed &
                       )

    ! update rho

    if ( fixed.eq.0 ) then

    dga = sum( ga(1:dmax) )
    call rnggamma( vrbA, dble(dga)+a_rho)
    call rnggamma( vrbB, dble(dmax-dga)+b_rho)
    rho =  vrbA/( vrbA +vrbB )

    end if

    ! update sig2

    dga = 0
    BtXtY = 0.0d0
    BtXtXB = 0.0d0
    sum_beta_sq = 0.0d0
    do j = 1, dmax
        if ( ga(j).eq.1 ) then
            ! update dga
            dga = dga +1
            ! update BtXtY
            BtXtY = BtXtY +beta(j)*XtY(j)
            ! update BtXtXB
            vrb = 0.0d0
            do i = 1, dmax
                if ( ga(i).eq.1 ) then
                    vrb = vrb + XtX(i,j)*beta(i)
                end if
            end do
            BtXtXB = BtXtXB + vrb*beta(j)
            ! update sum_beta_sq
            sum_beta_sq = sum_beta_sq +beta(j)**2
        end if
    end do

    vrbA = 0.5d0*y_en +0.5d0*dga +a_sig2
    vrbB = 0.5d0*( YtY +BtXtXB -2.d0*BtXtY ) &
            +lambda*sum_beta_sq &
            +b_sig2

    vrb = sig2
    call rnggamma( sig2, vrbA)
    sig2 = sig2/vrbB
    sig2 = 1.d0/sig2

    if ( sig2.ge.sig2_thr ) sig2 = vrb

    ! update lambda

    vrbA = 0.5d0*dble(dga)+a_lambda
    vrbB = 0.5d0*sum_beta_sq/sig2+b_lambda
    call rnggamma( lambda, vrbA)
    lambda = lambda/vrbB

end do

end subroutine mcmcsweep

! update ga and beta ==========================================================

subroutine update_ga_beta(Pr_ga,ga, &
                            rho, &
                            beta, &
                            sig2, &
                            lambda, &
                            dmax, XtY, XtX, &
                            fixed &
                           )

    implicit none

    integer, intent(in)             :: dmax
    integer, intent(inout)          :: ga(dmax)
    double precision, intent(in)    :: rho
    double precision, intent(inout) :: beta(dmax)
    double precision, intent(in)    :: sig2
    double precision, intent(in)    :: lambda
    double precision, intent(in)    :: XtY(dmax)
    double precision, intent(in)    :: XtX(dmax,dmax)
    integer, intent(in)             :: fixed
    double precision, intent(out)   :: Pr_ga(dmax)

!f2py intent(in)        :: dmax
!f2py intent(in,out)    :: ga
!f2py intent(in)        :: rho
!f2py intent(in,out)    :: beta
!f2py intent(in)        :: sig2
!f2py intent(in)        :: lambda
!f2py intent(in)        :: XtY
!f2py intent(in)        :: XtX
!f2py intent(in)        :: fixed
!f2py intent(out)       :: Pr_ga

    double precision, parameter     :: pi = 3.1415926535897932384626433832795d0

    integer                         :: j
    double precision                :: sj2
    double precision                :: mj
    double precision                :: rnd
    integer                         :: ivec(dmax)
    integer                         :: jj

    ! generate the permutation

    do jj = 1,dmax
        ivec(jj) = jj
    end do

    call genprm(ivec(1:dmax),dmax)

    do jj = 1, dmax

        ! pick a gPC bases -----------------------------------------------------

        j = ivec(jj)

        ! compute quantities --------------------------------------------------

        mj = ( &
                XtY(j) &
                -sum(XtX(1:dmax,j)*beta(1:dmax)) &
                +XtX(j,j)*beta(j) &
                ) / ( XtX(j,j) + lambda )

        sj2 = sig2 / ( XtX(j,j) + lambda )

        ! update ga -----------------------------------------------------------

        if ( fixed.eq.0 ) then

            if ( j .eq. 1  ) then

                ga(j) = 1
                Pr_ga(j) = 1.0

            else

                ! ... compute the probability

                Pr_ga(j) = 1.0d0/( &
                        1.0d0 &
                        +(1.d0-rho)/rho &
                            *sqrt(XtX(j,j)+lambda)/sqrt(lambda) &
                                *exp( -0.5*mj**2/sj2 ) )

                ! ... upgate ga
                call rnguniform( rnd )
                if ( Pr_ga(j).ge.rnd ) then
                    ga(j) = 1
                else
                    ga(j) = 0
                end if

            end if

        end if

        ! upgate beta ---------------------------------------------------------

        if (ga(j) .eq. 0) then
            beta(j) = 0.0d0
        else
            call rngnormal(beta(j:j),1)
            beta(j) = mj +beta(j)*sqrt(sj2)
        end if


!IF (SUM(GA).GT.10) EXIT

    end do

end subroutine update_ga_beta


