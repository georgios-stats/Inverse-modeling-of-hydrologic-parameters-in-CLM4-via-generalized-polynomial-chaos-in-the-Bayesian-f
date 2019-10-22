! ===============================
! GEORGIOS KARAGIANNIS
! GEORGIOS.KARAGIANNIS@PNNL.GOV
! PNNL
! 2013-08-20
! ===============================

subroutine Generate_sufficient_data(YtY,XtY,XtX,y_en,dmax,X,y)

	implicit none

	integer, intent(in)				:: y_en
	integer, intent(in)				:: dmax
	double precision, intent(in)	:: y(y_en)
	double precision, intent(in)	:: X(y_en,dmax)

	double precision, intent(out)	:: YtY
	double precision, intent(out)	:: XtY(dmax)
	double precision, intent(out)	:: XtX(dmax,dmax)

	integer								:: i
	integer								:: j
	integer								:: k
	double precision					:: XtX_ij
	double precision					:: XtY_j

	do j = 1, dmax
		do i = 1, j
			! XtX
			XtX_ij = 0.0d0
			do k = 1, y_en
				XtX_ij = XtX_ij + X(k,i)*X(k,j)
			end do
			XtX(i,j) = XtX_ij
            XtX(j,i) = XtX_ij
		end do 
		! XtY
		XtY_j = 0.0d0
		do k = 1, y_en
			XtY_j = XtY_j + X(k,j)*y(k)
		end do
		XtY(j) = XtY_j
	end do

	! YtY
	YtY = 0.0d0
	do k = 1, y_en
		YtY = YtY + y(k)*y(k)
	end do

end subroutine




