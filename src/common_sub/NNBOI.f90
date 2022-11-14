!-----------------------------------------------------------------------
!            This SUBROUTINE fill a matrix with $V_{n,n}$              !
!       at the Bohr-Oppenheimer Approximation, Coulomb particles       !
! na        :: number of atoms                                         !
! dis_mat   :: distance between nuclei matrix                          !
! a_num     :: atomi number (charge nucleus)                           !
! E_cmat    :: coulomb energy between nucleus in a matrix              !
! E_cou     :: sum coulomb energy at system                            !
!-----------------------------------------------------------------------
SUBROUTINE NNBOI(na, dis_mat, a_num, Ennmat, Enn)

IMPLICIT NONE
INTEGER                        :: i, j, na
INTEGER, DIMENSION(na)         :: a_num
REAL(KIND=8)                   :: Enn
REAL(KIND=8), EXTERNAL         :: f_coulomb
REAL(KIND=8), DIMENSION(na,na) :: dis_mat, Ennmat

Enn = 0.d0; Ennmat(:,:) = 0.d0

DO i = 1, na-1
  DO j = i+1, na
    Ennmat(i,j) = f_coulomb(a_num(i), a_num(j), dis_mat(i,j))
    Ennmat(j,i) = Ennmat(i,j)
    Enn         = Enn + Ennmat(i,j)
  ENDDO
ENDDO

ENDSUBROUTINE NNBOI
