!-----------------------------------------------------------------------
!   This SUBROUTINE gets the data about the basis set for the system   !
!   Right now just: 3GTO & 6-31G                                       !
!* basis_set :: basis set name, should be one of the already developd  *
!* na        :: nomber of atoms                                        *
!* atom      :: chemical symbolf of every atom                         *
!* bf        :: how many basis functions                               *
!* a_num     :: atomic number                                          *
!* n_bf_pa   :: number of basis functions per atom                     *
!* n_pri_bf  :: number of primitives per basis function                *
!* max_bf    :: maximum number of basis function                       *
!* max_prim  :: maximum number of primitives                           *
!* max_z     :: maximum number of zeta at basis for the system         *
!* zeta      :: $\zeta$ values for gaussian functions                  *
!* d         :: $d_{i,j}$ NORMALIZED                                   *
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BASIS_AT(basis_set, na, atom, bf, a_num, n_bf_pa, n_pri_bf, max_bf, max_prim, max_zeta, zeta, d)
IMPLICIT NONE
INTEGER                                     :: i, j, k, l, na, bf, info, max_bf, max_prim, max_zeta
INTEGER, DIMENSION(na)                      :: a_num, n_bf_pa
INTEGER, DIMENSION(max_zeta)                :: n_pri_bf
INTEGER, DIMENSION(:,:), ALLOCATABLE        :: n_prim
REAL(KIND=8), DIMENSION(max_zeta,max_prim)  :: zeta, d
REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: zeta_atom, d_atom
CHARACTER(LEN=15)                           :: basis_set
CHARACTER(LEN=2), DIMENSION(na)             :: atom
CHARACTER(LEN=30), DIMENSION(na)            :: file_basis_set
!-----------------------------------------------------------------------------------------------------------------------------------

!                                     (at atom, at basis, at Gaussian primitive)
ALLOCATE(n_prim(na, max_bf), zeta_atom(na, max_bf, max_prim), d_atom(na, max_bf, max_prim))

DO i = 1, na
  file_basis_set(i) = "../basis/"//TRIM(basis_set)//"/"//TRIM(atom(i))//".dat"
ENDDO

zeta(:,:) = 0.d0; d(:,:) = 0.d0
bf = 0; l = 1
DO i = 1, na
  OPEN(12,FILE=file_basis_set(i))
    READ(12,*) a_num(i), n_bf_pa(i)
!   (at the i'th atom, at the j'th basis of the i'th atom)
    DO j = 1, n_bf_pa(i)
      READ (12,*) n_prim(i,j)
      n_pri_bf(l)  = n_prim(i,j)
!     (at the i'th atom, at the j'th basis of the i'th atom, at th k'th zeta/d
      DO k = 1, n_prim(i,j)
        READ(12,*) zeta_atom(i,j,k), d_atom(i,j,k)
        zeta(l,k)  = zeta_atom(i,j,k)
        d(l,k)     = d_atom(i,j,k)
      ENDDO
      l = l + 1
    ENDDO
  CLOSE(12)
  bf = bf + n_bf_pa(i)
ENDDO
l = l - 1

IF (l .NE. bf) info = 211       !counters must be equal

CALL WIZARD(info)                                         !Everthing ok?

ENDSUBROUTINE BASIS_AT
