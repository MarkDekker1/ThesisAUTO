!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   Brusselator
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

	F(1) = PAR(1)-(PAR(2)+1)*U1+U1**2*U2
	F(2) = PAR(2)*U1-U1**2*U2

      END SUBROUTINE FUNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE STPNT(NDIM,U,PAR,Z)
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: Z

      DOUBLE PRECISION pi,A,B,Dx,Dy,RL

      pi=4*DATAN(1.d0)

!      *Set the (constant) parameters
       A  = 2.
       B  = 5.45
       Dx = 0.008
       Dy = 0.004
       RL = 0.75

       PAR(1)=A
       PAR(2)=B
       PAR(3)=RL

!      *Set the actual width of the space interval [0,PAR(11)]
       PAR(11) = 1.

!      *Set the initial data in the (scaled) interval [0,1]
       U(1) = A   - 0.5*DSIN(pi*Z)
       U(2) = B/A + 0.7*DSIN(pi*Z)

!      *Also set the space-derivative of the initial data
!      *Note the scaling by PAR(11)
       U(3) = - 0.5*pi*DCOS(pi*Z)/PAR(11)
       U(4) =   0.7*pi*DCOS(pi*Z)/PAR(11)

!      *Set the diffusion constants
       PAR(15) = Dx/RL**2
       PAR(16) = Dy/RL**2

      END SUBROUTINE STPNT

      SUBROUTINE BCND(NDIM,PAR,ICP,NBC,U0,U1,FB,IJAC,DBC) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), NBC, IJAC
      DOUBLE PRECISION, INTENT(IN) :: PAR(*), U0(NDIM), U1(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: FB(NBC)
      DOUBLE PRECISION, INTENT(INOUT) :: DBC(NBC,*)

!      *Define the boundary conditions (Dirichlet, in this example).
       FB(1)=U0(1)-PAR(1) 
       FB(2)=U0(2)-PAR(2)/PAR(1) 
       FB(3)=U1(1)-PAR(1) 
       FB(4)=U1(2)-PAR(2)/PAR(1) 

      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

      SUBROUTINE PVLS
      END SUBROUTINE PVLS
