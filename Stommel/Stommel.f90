!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   ab :            The A --> B reaction 
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

! Evaluates the algebraic equations or ODE right hand side

! Input arguments :
!      NDIM   :   Dimension of the ODE system 
!      U      :   State variables
!      ICP    :   Array indicating the free parameter(s)
!      PAR    :   Equation parameters

! Values to be returned :
!      F      :   ODE right hand side values

! Normally unused Jacobian arguments : IJAC, DFDU, DFDP (see manual)

      IMPLICIT NONE
      INTEGER NDIM, IJAC, ICP(*)
      DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)
      DOUBLE PRECISION U1,U2,eta1, eta2, eta3

	eta1 = PAR(1)
	eta2 = PAR(2)
	eta3 = PAR(3)

	U1 = U(1)
	U2 = U(2)

	F(1) = eta1 - U1*(1 + ABS(U1 - U2))
	F(2) = eta2 - U2*(eta3 + ABS(U1 - U2))

      END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

! Input arguments :
!      NDIM   :   Dimension of the ODE system 

! Values to be returned :
!      U      :   A starting solution vector
!      PAR    :   The corresponding equation-parameter values
!      T      :	  Not used here

      IMPLICIT NONE
      INTEGER NDIM
      DOUBLE PRECISION U(NDIM), PAR(*), T

! Initialize the equation parameters
       PAR(1)=2.
       PAR(2)=0.
       PAR(3)=0.3

! Initialize the solution
       U(1)=1.0
       U(2)=0.0

      END SUBROUTINE STPNT
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! The following subroutines are not used here,
! but they must be supplied as dummy routines

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
