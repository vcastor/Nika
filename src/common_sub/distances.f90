!-----------------------------------------------------------------------
! This SUBROUTINE fill matrices with the atom and basis function       !
! distances                                                            !
!* na      :: number of atoms                                          *
!* bf      :: number of basis functions                                *
!* nbfpa   :: number of basis functions per atom                       *
!* xyz     :: coordinates of atoms (bohr)                              *
!* xyz_bf  :: coordinates of basis function (bohr)                     *
!* dis_a   :: distance between atoms                                   *
!* dis2_a  :: squared distance between atoms                           *
!* dis_bf  :: distance between basis functions                         *
!* dis2_bf :: squareddistance between basis functions                  *
!-----------------------------------------------------------------------
SUBROUTINE DISTANCES(na, bf, nbf_pa, xyz, xyz_bf, dis_a, dis2_a, &
                                                        dis_bf, dis2_bf)
IMPLICIT NONE
INTEGER                        :: i, j, k, l, m, n, na, bf
INTEGER, DIMENSION(na)         :: nbf_pa
REAL(KIND=8), EXTERNAL         :: f_distance2
REAL(KIND=8), DIMENSION(3,na)  :: xyz
REAL(KIND=8), DIMENSION(3,bf)  :: xyz_bf
REAL(KIND=8), DIMENSION(na,na) :: dis_a, dis2_a
REAL(KIND=8), DIMENSION(bf,bf) :: dis_bf, dis2_bf

DO i = 1, na                                     !distance between atoms
  DO j = i, na
    dis2_a(i,j) = f_distance2(xyz(:,i), xyz(:,j))
    dis2_a(j,i) = dis2_a(i,j)
    dis_a(i,j)  = DSQRT(dis2_a(i,j))
    dis_a(j,i)  = dis_a(i,j)       !same distance between A & B as B & A
  ENDDO
ENDDO

k = 1                                                           !counter
DO i = 1, na                            !coordinates for basis functions
  DO j = 1, nbf_pa(i)
    xyz_bf(:,k) = xyz(:,i)
    k = k + 1
  ENDDO
ENDDO

n = 1; m = 1                                                   !counters
DO i = 1, na                          !distances between basis functions
  DO j = 1, na
    DO k = 1, nbf_pa(i)
      DO l = 1, nbf_pa(j)
        dis2_bf(n,m) = dis2_a(j,i)
        dis_bf(n,m)  = dis_a(j,i)
        n = n + 1                                     !only nbf_pa times
      ENDDO
      n = n - nbf_pa(j)     !restart every time that we jump to next row
      m = m + 1
    ENDDO
    n = n + nbf_pa(j)   !start where we stop last time with the new vale
    m = m - nbf_pa(i)    !restart every time that we jump to next column
  ENDDO
  n = 1                                                       !every row
  m = m + nbf_pa(i)
ENDDO

ENDSUBROUTINE DISTANCES
