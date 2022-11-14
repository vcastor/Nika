!-----------------------------------------------------------------------------------------------------------------------------------
!********                                    UNRESTRICTED  HARTREE-FOCK-ROOTHAN ALRGORITHM                                  ********
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE UHF(bf, n_occ, S, T, V, TwoEleInt, E_ele)
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER                               :: bf, n_occ, na, nb, i, j, k, l, m, n, step
REAL (KIND=8)                         :: threshold, aux1, aux2, E_ele
REAL (KIND=8), DIMENSION(bf)          :: ea, eb
REAL (KIND=8), DIMENSION(bf,bf)       :: S, T, V, P, Pa, Pb, P_aux, Fa, Fb, Hcore, Jm, Kam, Kbm, Fap, Fbp, Ca, Cb, Cap, Cbp, X
REAL (KIND=8), DIMENSION(bf,bf,bf,bf) :: TwoEleInt
!-----------------------------------------------------------------------------------------------------------------------------------

na = n_occ 
nb = n_occ + 1

Hcore(:,:) = T(:,:) + V(:,:)                              !core Hamiltonian: H^{core} = KineticIntegral + Nuclear Atraction Integral
Ca(:,:) = 0.d0; Cb(:,:) = 0.d0; P_aux(:,:)= 1.d0;                                      !just in case the RAM was with something else

!-----------------------------------------------------------------------------------------------------------------------------------
!---- Here start the loop until the threshold smaller than $10^{-10}$

threshold = 1.d0; step = 0
DO WHILE (threshold .GT. 1d-10)

    Pa(:,:) = 0.d0; Pb(:,:) = 0.d0       !Restart the sum every loop
    DO i = 1, bf                             !Density matrix: P(i,j)
      DO j = 1, bf
        DO k = 1, na
          Pa(j,i) = Pa(j,i) + Ca(j,k)*Ca(i,k)
        ENDDO
        DO k = 1, nb
          Pb(j,i) = Pb(j,i) + Cb(j,k)*Cb(i,k)
        ENDDO
      ENDDO
    ENDDO

    P = Pa + Pb                      !P^{Total} = P^{\alpha} + P^{\beta}

    Jm(:,:) = 0.d0; Kam(:,:) = 0.d0; Kbm(:,:) = 0.d0 !restart every loop
    DO n = 1, bf                         !J = P * (\mu\nu|\lambda\sigma)
      DO m = 1, bf                       !K = P * (\mu\lambda|\nu\sigma)
        DO l = 1, bf
          DO k = 1, bf                   !k :: \sigma
            Jm(m,n)  = Jm(m,n)  +  P(k,l)*TwoEleInt(m,n,l,k)
            Kam(m,n) = Kam(m,n) + Pa(k,l)*TwoEleInt(m,l,n,k)
            Kbm(m,n) = Kbm(m,n) + Pb(k,l)*TwoEleInt(m,l,n,k)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

!   F^\alpha = T + V + P^T(\mu\nu|\lambda\sigma) - P^\alpha(\mu\lambda|\nu\sigma)
!   F^\beta  = T + V + P^T(\mu\nu|\lambda\sigma) - P^\beta(\mu\lambda|\nu\sigma)
    Fa = Hcore + Jm - Kam
    Fb = Hcore + Jm - Kbm

!***********************************************************************
!********                   F^αC^α = SC^α ε^α                   ********
!********                   F^βC^β = SC^β ε^β                   ********
!********                         ⬇️                             ********
!                                                                      !
    CALL ORTHOGONALIZATION_M(bf, S, X)           !      X = S^(-1/2)   !
    Fap = MATMUL(TRANSPOSE(X),MATMUL(Fa,X))      !    F'α = X^T Fα X   !
    Fbp = MATMUL(TRANSPOSE(X),MATMUL(Fb,X))      !    F'β = X^T Fβ X   !
    Cap = Fap; CALL DIAG_M(bf, Cap, ea)          ! compute C'α and εα  !
    Cbp = Fbp; CALL DIAG_M(bf, Cbp, eb)          ! compute C'β and εβ  !
    Ca  = MATMUL(X,Cap)                          !      C'α ➡️   Cα     !
    Cb  = MATMUL(X,Cbp)                          !      C'β ➡️   Cβ     !
!                                                                      !
!********                         ⬇️                             ********
!********                 F'^βC'^β = C'^βε^β                    ********
!********                 F^α'C'^α = C'^αε^α                    ********
!***********************************************************************

    aux1 = 0.d0; aux2 = 0.d0                     !            Converged?
    DO j = 1, bf                                 
      DO i = 1, bf
        aux1 = P(i,j) - P_aux(i,j)
        aux2 = aux2 + aux1*aux1
      ENDDO
    ENDDO
    threshold = DSQRT(aux2/REAL(bf*bf))

    P_aux = P                         !Save the Density Matrix of the steps

    step = step + 1                                        !Step counter
ENDDO                                                                   !Close the DO WHILE of the Hartree-Fock iterative part 

!-----------------------------------------------------------------------------------------------------------------------------------
!---- Finally. Compute the energy at the converged system

CALL EUHF(bf, P, T, V, Jm, Kam, Kbm, E_ele)

ENDSUBROUTINE
!-----------------------------------------------------------------------------------------------------------------------------------
!********                                       END OF THE HARTREE-FOCK-ROOTHAN ALGORITHM                                   ********
!-----------------------------------------------------------------------------------------------------------------------------------
