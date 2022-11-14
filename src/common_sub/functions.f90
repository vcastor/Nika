!***********************************************************************
!********************   HONI SOIT QUI MAL Y PENSE   ********************
!********************        NOLI SE TANGERE        ********************
!***********************************************************************
    FUNCTION f_distance2(xyz, uvw) RESULT (distance2)
    IMPLICIT NONE
    REAL (KIND=8)               :: distance2, dx, dy, dz
    REAL (KIND=8), DIMENSION(3) :: xyz, uvw

    dx = xyz(1) - uvw(1)
    dy = xyz(2) - uvw(2)
    dz = xyz(3) - uvw(3)

    distance2 = dx*dx + dy*dy + dz*dz

    ENDFUNCTION f_distance2
!***********************************************************************
    FUNCTION f_distance(xyz, uvw) RESULT (distance)
    IMPLICIT NONE
    REAL (KIND=8)               :: distance
    REAL (KIND=8), EXTERNAL     :: f_distance2
    REAL (KIND=8), DIMENSION(3) :: xyz, uvw

    distance = DSQRT(f_distance2(xyz, uvw))

    ENDFUNCTION f_distance
!***********************************************************************
    FUNCTION f_coulomb(q,Z,distance) RESULT (coulomb)
    INTEGER       :: q, Z
    REAL (KIND=8) :: coulomb, distance

    coulomb = DBLE(q*Z)/distance

    ENDFUNCTION f_coulomb
!***********************************************************************
    FUNCTION f_cero(x) RESULT (y)
    IMPLICIT NONE
    INCLUDE 'parameters.h'
    REAL (KIND=8), INTENT(IN) :: x
    REAL (KIND=8)             :: y

    IF (x .GT. 1.d-10) THEN                                 !Aproximation
      y  = 0.5d0*DSQRT(pi/x)*DERF(DSQRT(x))
    ELSE                                                    !Analytically
      y = 1.d0
    ENDIF

    ENDFUNCTION f_cero
!***********************************************************************
    FUNCTION f_trace_mat(n, A) RESULT (tr)
    IMPLICIT NONE
    INTEGER                       :: i, n
    REAL (KIND=8)                 :: tr
    REAL (KIND=8), DIMENSION(n,n) :: A

    tr = 0.d0
    DO i = 1, n
      tr = tr + A(i,i)
    ENDDO

    ENDFUNCTION
!***********************************************************************
