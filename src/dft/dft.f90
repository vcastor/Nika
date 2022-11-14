!-----------------------------------------------------------------------------------------------------------------------------------
!********                                 RESTRICTED  DENSITY FUNCTIONAL THEORY ALRGORITHM                                  ********
!-----------------------------------------------------------------------------------------------------------------------------------
! bf        :: number of basis functions                                                                                           !
! n_occ     :: number of occupied orbitals                                                                                         !
! S         :: Overlap integral                                                                                                    !
! T         :: Kinetic integral                                                                                                    !
! V         :: Potential integral                                                                                                  !
! TwoEleInt :: Two electron integral                                                                                               !
! exchf     :: exchange functional                                                                                                 !
! corrf     :: correlation functionla                                                                                              !
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE DFT(bf, n_occ, S, T, V, TwoEleInt, exchf, corrf, E_ele)
IMPLICIT NONE
INTEGER                              :: bf, n_occ
REAL(KIND=8)                         :: E_ele
REAL(KIND=8), DIMENSION()
REAL(KIND=8), DIMENSION(bf,bf)       :: S, T, V
REAL(KIND=8), DIMENSION(bf,bf,bf,bf) :: TwoEleInt
CHARACTER(LEN=)                     :: exchf, corrf
ENDSUBROUTINE DFT
