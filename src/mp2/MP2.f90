!-----------------------------------------------------------------------------------------------------------------------------------
!********                                             RESTRICTED  MP2  ALRGORITHM                                           ********
!-----------------------------------------------------------------------------------------------------------------------------------
! bf        :: number of basis functions                                                                                           !
! n_occ     :: number of occupied orbitals                                                                                         !
! S         :: Overlap integral                                                                                                    !
! T         :: Kinetic integral                                                                                                    !
! V         :: Potential integral                                                                                                  !
! TwoEleInt :: Two electron integral                                                                                               !
! E_ele     :: Electronic energy (Hartree-Fock)                                                                                    !
! E_cmp1    :: Electronic correlation energy (MP2)                                                                                 !
! Cmat      :: Coefficent matrix, from Eigenvalue problem (FCS=CE)                                                                 !
! AO        :: Atomic Orbital two electron integral (as tensor)                                                                    !
! MO        :: Molecular Orbital two electron integral (as tensor)                                                                 !
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE MP2(bf, n_occ, S, T, V, TwoEleInt, E_hf, E_cmp2, step)
!-----------------------------------------------------------------------------------------------------------------------------------
IMPLICIT NONE
INTEGER                                :: bf, n_occ, i, j, k, l
REAL (KIND=8)                          :: E_hf, E_cmp2, eps, E_ele, step
REAL (KIND=8), DIMENSION(bf)           :: e
REAL (KIND=8), DIMENSION(n_occ)        :: h
REAL (KIND=8), DIMENSION(bf,bf)        :: S, T, V, Hcore, C
REAL (KIND=8), DIMENSION(bf,bf,bf,bf)  :: TwoEleInt, MO

!-----------------------------------------------------------------------
!--- First call Hartree-Fock-Roothaan Algorithm
CALL RHF(bf, n_occ, S, T, V, TwoEleInt, E_ele, C, e, step)

!-----------------------------------------------------------------------
!--- We need the two electron integrals in Molecular Orbitals
CALL AOtoMO(bf,C,TwoEleInt,MO)

!-----------------------------------------------------------------------
!--- Compute with the Core Hamiltonian
Hcore = T + V
h(:)  = 0.d0 
DO i = 1, n_occ
  DO j = 1, bf
    DO k = 1, bf
      !               μ          μ ν    ν
      h(i) = h(i) + C(j,i)*Hcore(j,k)*C(k,i)
    ENDDO
  ENDDO
ENDDO

!-----------------------------------------------------------------------
!--- Hartree Fock Energy with the Molecular Orbitals
E_hf = 0.d0
DO i = 1, n_occ
  E_hf = E_hf + 2.d0*h(i)
  DO j = 1, n_occ
    E_hf = E_hf + 2.d0*MO(i,j,i,j) - MO(i,j,j,i)
  ENDDO
ENDDO

!-----------------------------------------------------------------------
!---                     MP2 correlation energy                      ---
E_cmp2 = 0.d0
DO i = 1, n_occ
  DO j = 1, n_occ
    DO k  = n_occ+1, bf
      DO l = n_occ+1, bf
        eps    = e(i) + e(j) - e(k) - e(l)
        E_cmp2 = E_cmp2 + MO(i,j,k,l)*(2.d0*MO(i,j,k,l) - MO(i,j,l,k))/eps
      ENDDO
    ENDDO
  ENDDO
ENDDO

ENDSUBROUTINE MP2
