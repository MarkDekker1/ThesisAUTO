!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   Rombouts Vegetation model
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----


      IMPLICIT NONE
      INTEGER NDIM, IJAC, ICP(*)
      DOUBLE PRECISION U(NDIM), PAR(*), F(NDIM), DFDU(*), DFDP(*)
      DOUBLE PRECISION U1,U2

	U1 = U(1)
	U2 = U(2)

	F(1) = PAR(2)*U1*(1-U1)-U1*U2-PAR(1)*(1-EXP(-PAR(3)*U1))
	F(2) = -U2+PAR(4)*U1*U2

      END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)

      IMPLICIT NONE
      INTEGER NDIM
      DOUBLE PRECISION U(NDIM), PAR(*), T

! Initialize the equation parameters
       PAR(1)=0.	!a
       PAR(2)=3.	!b
       PAR(3)=5.	!c
       PAR(4)=3.	!d

! Initialize the solution
       U(1)= 0.
       U(2)= 0.

      END SUBROUTINE STPNT
!----------------------------------------------------------------------
!----------------------------------------------------------------------
      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

SUBROUTINE PVLS(NDIM,U,PAR)
!--------- ----

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NDIM
  DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
  DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
  DOUBLE PRECISION GETP

! Set PAR(9) equal to U1 
  PAR(9)=GETP('BV0',1,U)

END SUBROUTINE PVLS
!----------------------------------------------------------------------
!----------------------------------------------------------------------
