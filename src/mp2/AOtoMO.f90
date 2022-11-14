!-----------------------------------------------------------------------
! This SUBROUTINE computes:                                            !
!                           (μν|λσ) ➡️  (pq|rs)                         !
!         Changing from Atomic Orbitals to Molecular Orbitals          !
! bf   :: how many basis functions                                     !
! Cmat :: Coefficent matrix, from Eigenvalue problem (FCS=CE)          !
! AO   :: Atomic Orbital two electron integral (as tensor)             !
! MO   :: Molecular Orbital two electron integral (as tensor)          !
!-----------------------------------------------------------------------

SUBROUTINE AOtoMO(bf,Cmat,AO,MO)

IMPLICIT NONE
INTEGER                              :: i, m, n, l, q, r, s, bf
REAL(KIND=8), DIMENSION(bf,bf)       :: Cmat
REAL(KIND=8), DIMENSION(bf,bf,bf,bf) :: AO, MO, tmp

! 1st step (but with the 4th term)
tmp(:,:,:,:) = 0.d0
DO s = 1, bf
  DO i = 1, bf
    DO l = 1, bf
      DO n = 1, bf
        DO m = 1, bf
!            (μ ν λ s)      (μ ν λ s)       (σ s)   (μ ν λ σ)
          tmp(m,n,l,i) = tmp(m,n,l,i) + Cmat(s,i)*AO(m,n,l,s)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO

! 2nd step (with 3rd term)
MO(:,:,:,:) = 0.d0
DO l = 1, bf
  DO s = 1, bf
    DO i = 1, bf
      DO n = 1, bf
        DO m = 1, bf
!           (μ ν r s)     (μ ν r s)       (λ r)    (μ ν λ s)
          MO(m,n,i,s) = MO(m,n,i,s) + Cmat(l,i)*tmp(m,n,l,s)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO

! 3rd step (with 2nd term)
tmp(:,:,:,:) = 0.d0
DO n = 1, bf
  DO s = 1, bf
    DO r = 1, bf
      DO i = 1, bf
        DO m = 1, bf
!            (μ q r s)      (μ q r s)       (ν q)   (μ ν r s)
          tmp(m,i,r,s) = tmp(m,i,r,s) + Cmat(n,i)*MO(m,n,r,s)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO

! 4th step (with 1st term)
MO(:,:,:,:) = 0.d0
DO m = 1, bf
  DO s = 1, bf
    DO r = 1, bf
      DO q = 1, bf
        DO i = 1, bf
!           (p q r s)     (p q r s)       (μ p)    (μ q r s)
          MO(i,q,r,s) = MO(i,q,r,s) + Cmat(m,i)*tmp(m,q,r,s)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO

ENDSUBROUTINE AOtoMO
