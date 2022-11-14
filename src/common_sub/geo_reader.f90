!-----------------------------------------------------------------------
! This SUBROUTINE reads the coordiantes in:                            !
!             input file for Nika proyect                              !
! input_name :: file input name                                        !
! basis_set  :: basis set name                                         !
! na         :: number of atoms                                        !
! atom       :: atom chemical characer                                 !
! xyz        :: coordinates                                            !
! info       :: info flag                                              !
!-----------------------------------------------------------------------
SUBROUTINE GEO_READER(input_name, na, atom, xyz, units)

IMPLICIT NONE
INCLUDE 'parameters.h'
INTEGER                         :: i, j, na, stat, ic, info
REAL(KIND=8), DIMENSION(3,na)   :: xyz
CHARACTER(LEN=15)               :: units
CHARACTER(LEN=100)              :: input_name, buf
CHARACTER(LEN=2), DIMENSION(na) :: atom

OPEN(11, FILE=input_name, ACTION='READ')
  ic = 0
  DO
    READ(11, *, IOSTAT=stat) buf
    IF (stat.NE.0) EXIT
    ic = ic + 1
  ENDDO
CLOSE(11)

IF (na.EQ.(ic-6)) info = 0
IF (na.GT.(ic-6)) info = 104
IF (na.LT.(ic-6)) info = 105
CALL WIZARD(info)                                         !Everthing ok?

OPEN(12,FILE=input_name)

  DO i = 1, 6; READ(12,*); ENDDO         !dummy lines

  DO i = 1, na                           !chemical symbol, coordinates
    READ(12,*) atom(i), (xyz(j,i), j=1, 3)
  ENDDO

CLOSE(12)

SELECTCASE (units)
  CASE ('angstrom')
    xyz(:,:) = atb*xyz(:,:)
ENDSELECT

ENDSUBROUTINE GEO_READER
