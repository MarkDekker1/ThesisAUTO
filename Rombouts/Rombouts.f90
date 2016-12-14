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
      DOUBLE PRECISION U1,U2, beta,alpha,R0,alpha0

	U1 = U(1)
	U2 = U(2)

	alpha0 	= MAX(PAR(6)+MIN(0.,(PAR(7)-PAR(6))*(U1-PAR(8))/(PAR(9)-PAR(8))),PAR(7))
	alpha 	= (1-PAR(3))*alpha0+PAR(3)*(PAR(4)*U2+PAR(5)*(1-U2))
	beta	= MAX(0.,1-PAR(13)*(U1-PAR(12))**2)
	R0	= PAR(10) + PAR(11)*(U1 - PAR(12))

	F(1) = ((1 - alpha)*PAR(2) - R0)/PAR(1)
	F(2) = beta*U2*(1 - U2) - PAR(14)*U2

      END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,T)

      IMPLICIT NONE
      INTEGER NDIM
      DOUBLE PRECISION U(NDIM), PAR(*), T

! Initialize the equation parameters
       PAR(1)=500.0	!Ct
       PAR(2)=342.5	!Q0
       PAR(3)=0.3	!p
       PAR(4)=0.1	!av
       PAR(5)=0.4	!ag
       PAR(6)=0.85	!alphamax
       PAR(7)=0.25	!alphamin
       PAR(8)=263.0	!Tal
       PAR(9)=300.0	!Tau
       PAR(10)=200.0	!B0
       PAR(11)=2.5	!B1
       PAR(12)=283.0	!Topt
       PAR(13)=0.004	!k
       PAR(14)=3.4608467752E-01!0.41777072213	!gmma

! Initialize the solution
       U(1)= 2.9432679889E+02!290.761399859
       U(2)= 2.8908303020E-01!0.449608397214

      END SUBROUTINE STPNT
!----------------------------------------------------------------------
!----------------------------------------------------------------------
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
