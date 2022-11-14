!-----------------------------------------------------------------------
!   This SUBROUTINE write a file if we have troubles 
!***********************************************************************
!       with hate,
!               Victoria Castor 2021
!-----------------------------------------------------------------------

SUBROUTINE writer_error(info)
INTEGER :: info

OPEN(666,FILE="hfmp2dft.err")

!-----------------------------------------------------------------------------------------------------------------------------------
SELECTCASE (info)
  CASE(101)
    WRITE(666,*) "The only 2 units that the program can read are:"
    WRITE(666,*) "i)  angstorom"
    WRITE(666,*) "ii) bohr"
    WRITE(666,*) "working in more, like Imperial Units... "
    WRITE(666,*) "                                      HAHAHA! of course not"
    WRITE(666,*) "just that two above 😅"
    WRITE(666,*) '                                         SORRY NOT SORRY 😬'

  CASE(201)
    WRITE(666,*) 'There are someting wrong in the basis set data, PLEASE:'
    WRITE(666,*) ' -> contact to Victoria Castor [at] :'
    WRITE(666,*) ' '
    WRITE(666,*) 'e-mail:      vcastorv@gmail.com'
    WRITE(666,*) 'Instagram:   @vcastorv'
    WRITE(666,*) 'Twitter:     @vcastor'
    WRITE(666,*) '                                         SORRY NOT SORRY 😬'

  CASE(501)
    WRITE(666,*) 'Overlap integral numerical error'
    WRITE(666,*) 'look at: ./tmp/Overlap.int'

  CASE(502)
    WRITE(666,*) 'Potential integral numerical error'
    WRITE(666,*) 'look at: ./tmp/Potential.int'

  CASE(503)
    WRITE(666,*) 'Kinetic integral numerical error'
    WRITE(666,*) 'look at: ./tmp/Kinetic.int'

  CASE(504)
    WRITE(666,*) 'Two electron integral numerical error'
    WRITE(666,*) 'look at: ./tmp/TwoElectron.int'

  CASE(601)
    WRITE(666,*) 'Rotation numerical error'
    WRITE(666,*) 'Values too small to compute a Lowdin orthonormalization'

ENDSELECT
CLOSE(666)
ENDSUBROUTINE writer_error
