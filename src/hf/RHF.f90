!-----------------------------------------------------------------------------------------------------------------------------------
!********                                    RESTRICTED  HARTREE-FOCK-ROOTHAN ALRGORITHM                                    ********
!-----------------------------------------------------------------------------------------------------------------------------------
! bf        :: number of basis functions                                                                                           !
! n_occ     :: number of occupied orbitals                                                                                         !
! S         :: Overlap integral                                                                                                    !
! T         :: Kinetic integral                                                                                                    !
! V         :: Potential integral                                                                                                  !
! TwoEleInt :: Two electron integral                                                                                               !
! E_ele     :: Electronic energy                                                                                                   !
! C         :: C matrix LCAO                                                                                                       !
! e         :: e vector                                                                                                            !
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE RHF(bf, n_occ, S, T, V, TwoEleInt, E_ele, C, e, step)
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER                               :: bf, n_occ, i, j, k, l, m, n, step
REAL (KIND=8)                         :: threshold, aux1, aux2, E_ele
REAL (KIND=8), DIMENSION(bf)          :: e
REAL (KIND=8), DIMENSION(bf,bf)       :: S, T, V, C, P, Pb, F, Hcore, Jm, Km, Fp, Cp, X
REAL (KIND=8), DIMENSION(bf,bf,bf,bf) :: TwoEleInt
!-----------------------------------------------------------------------------------------------------------------------------------

Hcore(:,:) = T(:,:) + V(:,:)                              !core Hamiltonian: H^{core} = KineticIntegral + Nuclear Atraction Integral
C(:,:) = 0.d0; Pb(:,:)= 1.d0                                                           !just in case the RAM was with something else

!-----------------------------------------------------------------------------------------------------------------------------------
!---- Here start the loop until the threshold smaller than $10^{-10}$

threshold = 1.d0; step = 0
DO WHILE (threshold .GT. 1d-10)

    P(:,:) = 0.d0                            !Restart the sum every loop
    DO i = 1, bf                             !Density matrix: P(i,j)
      DO j = 1, bf
        DO k = 1, n_occ
          P(j,i) = P(j,i) + C(j,k)*C(i,k)
        ENDDO
      ENDDO
    ENDDO

    Jm(:,:) = 0.d0; Km(:,:) = 0.d0                   !restart every loop
    DO n = 1, bf                         !J = P * (\mu\nu|\lambda\sigma)
      DO m = 1, bf                       !K = P * (\mu\lambda|\nu\sigma)
        DO l = 1, bf
          DO k = 1, bf                   !k :: \sigma
            Jm(m,n) = Jm(m,n) + P(k,l)*TwoEleInt(m,n,l,k)
            Km(m,n) = Km(m,n) + P(k,l)*TwoEleInt(m,l,n,k)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

!   F = T + V + P(2(\mu\nu|\lambda\sigma)-(\mu\lambda|\nu\sigma))
    F = Hcore + 2.d0*Jm - Km                                            !F = H^{core} + 2J - K 

!***********************************************************************
!********                       FC = SCE                        ********
!********                         ⬇️                             ********
!                                                                      !
    CALL ORTHOGONALIZATION_M(bf, S, X)             !     X = S^(-1/2)  !
    Fp = MATMUL(TRANSPOSE(X),MATMUL(F,X))          !    F' = X^T F X   !
    Cp = Fp; CALL DIAG_M(bf, Cp, e)                ! compute C' and ε  !
    C  = MATMUL(X,Cp)                              !      C' ➡️   C     !
!                                                                      !
!********                         ⬇️                             ********
!********                     F'C' = C'ε                        ********
!***********************************************************************

    aux1 = 0.d0; aux2 = 0.d0                           !      Converged?
    DO j = 1, bf                                 
      DO i = 1, bf
        aux1 = P(i,j) - Pb(i,j)
        aux2 = aux2 + aux1*aux1
      ENDDO
    ENDDO
    threshold = DSQRT(aux2/REAL(bf*bf))

    Pb = P                         !Save the Density Matrix of the steps

    step = step + 1                                        !Step counter
ENDDO                                                                   !Close the DO WHILE of the Hartree-Fock iterative part 

!-----------------------------------------------------------------------------------------------------------------------------------
!---- Finally. Compute the energy at the converged system

CALL ERHF(bf, P, T, V, Jm, Km, E_ele)

ENDSUBROUTINE
!-----------------------------------------------------------------------------------------------------------------------------------
!********                                       END OF THE HARTREE-FOCK-ROOTHAN ALGORITHM                                   ********
!-----------------------------------------------------------------------------------------------------------------------------------
