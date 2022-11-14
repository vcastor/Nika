!-----------------------------------------------------------------------
!          This SUBROUTINE reads the input file for NIKA               !
!*                                               Works for all methods *
!* method    :: theory level: HF, MP2, DFT (keywords: hf, mp2, dft)    *
!* inp_n     :: input name, can be in other directory                  *
!* basis_set :: basis set name,                                        *
!* q         :: over charge of molecule                                *
!* na        :: number of atoms                                        *
!* info      :: info flag if something is wrong                        *
!-----------------------------------------------------------------------
SUBROUTINE OPT_READER(inp_n, method, basis_set, q, na, exchf, corrf)

IMPLICIT NONE
INCLUDE 'parameters.h'
INTEGER            :: q, na, info
CHARACTER(LEN=15)  :: method, basis_set, units, exchf, corrf
CHARACTER(LEN=100) :: inp_n

OPEN(11,FILE=inp_n)

    READ(11,*) method                                      !thoery level
    READ(11,*) basis_set                                 !basis set name
    READ(11,*) q                                            !over charge
    READ(11,*) units                                     !distance units
    READ(11,*) na                                       !number of atoms

CLOSE(11)

SELECTCASE (method)
    CASE('RHF')
        info = 0
    CASE('UHF')
        info = 0
    CASE('MP2')
        info = 0
    CASE('HF')
        method = 'UHF'
        info    = 0
    CASE('rhf')
        method = 'RHF'
        info   = 0
    CASE('uhf')
        method = 'UHF'
        info   = 0
    CASE('hf')
        method = 'UHF'
        info    = 0
    CASE('mp2')
        method = 'MP2'
        info   = 0
    CASE('BLYP')
        method = 'DFT'
        exchf  = 'Becke'
        corrf  = 'LYP'
        info   = 0
    CASE('blyp')
        method = 'DFT'
        exchf  = 'Becke'
        corrf  = 'LYP'
        info   = 0
    CASE('LDA')
        method = 'DFT'
        exchf  = 'Slater'
        corrf  = 'Wigner'
        info   = 0
    CASE('lda')
        method = 'DFT'
        exchf  = 'Slater'
        corrf  = 'Wigner'
        info   = 0
    !CASE()
    CASE DEFAULT
        info = 101
ENDSELECT

SELECTCASE (units)
    CASE ('angstrom')
        info = 0
    CASE ('bohr')
        info = 0
    CASE('Angstrom')
        units = 'angstrom'
        info  = 0
    CASE('Bohr')
        units = 'bohr'
        info  = 0
    CASE('ANGSTROM')
        units = 'angstrom'
        info  = 0
    CASE('BOHR')
        units = 'bohr'
        info  = 0
    CASE DEFAULT
        info = 102
ENDSELECT

IF (na .LE. 0) info = 103

CALL WIZARD(info)                                         !Everthing ok?

ENDSUBROUTINE OPT_READER
