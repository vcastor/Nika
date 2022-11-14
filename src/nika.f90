                            PROGRAM Nika
!***********************************************************************
! Nika proyect:      NIKA: Nika Isn't a Known Acronym                  !
!       -> Hartree-Fock-Roothaan                       (HFR)           !
!       -> Møller–Plesset at second order              (MP2)           !
!       -> Density Functional Thoery                   (DFT)           !
!***********************************************************************
! PLEASE DO NOT ATTEMP TO SIMPLIFY THIS CODE.                          !
! KEEP THE SPACE SHUTTLE FLYING.                                       !
!                     ~Just kidding, all contributions are wellcome <3 !
!***********************************************************************
! We call this style 'space shuttle style'. Space shuttle style is     !
! meant to ensure that every branch and condition is considered and    !
! accounted for the same way code is written at NASA for apllications  !
! like the space shuttle.                                              !
!***********************************************************************
! PLEASE: FIRST, READ THE README FILE                                  !
!***********************************************************************
! Full program was developed thinking in Hartree Atomic Units:         !
! Reduced Planck constant: ℏ = 1, (\hbar),       atomic unit of action !
! Elementary charge:       e = 1,                atomic unit of charge !
! Bohr radius:           a_0 = 1,                atomic unit of length !
! Electron mass:         m_e = 1                 atomic unit of mass   !
!***********************************************************************
!       with love,                                                     !
!            Victoria Castor 2022                                      !
!***********************************************************************
IMPLICIT NONE
!-----------------------------------------------------------------------------------------------------------------------------------
LOGICAL                                          :: closed_shell

INTEGER                                          :: i, q, na, bf, n_occ, n_vir, n_ele, max_prim, max_zeta, max_bf, step
INTEGER, DIMENSION(:), ALLOCATABLE               :: a_num, n_pri_bf, n_bf_pa, n_pri_bf_aux

REAL (KIND=8)                                    :: Enn, E_ele, E_cmp2, threshold
REAL (KIND=8), DIMENSION (:), ALLOCATABLE        :: e
REAL (KIND=8), DIMENSION (:,:), ALLOCATABLE      :: S, T, V, C, xyz, xyz_bf, zeta, d, Ennmat, dis_a, dis_bf, dis2_a, dis2_bf
REAL (KIND=8), DIMENSION (:,:,:,:), ALLOCATABLE  :: TwoEleInt

CHARACTER (LEN=2), DIMENSION (:), ALLOCATABLE    :: atom
CHARACTER (LEN=15)                               :: basis_set, units, method, exchf, corrf
CHARACTER (LEN=100)                              :: input_name
!-----------------------------------------------------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
!********                      INPUT DATA                       ********
!-----------------------------------------------------------------------

CALL GETARG(1, input_name)

CALL OPT_READER(input_name, method, basis_set, q, na, exchf, corrf)

ALLOCATE(atom(na), xyz(3,na))

CALL GEO_READER(input_name, na, atom, xyz, units)

!-----------------------------------------------------------------------
!********                      BASIS SET                        ********
!-----------------------------------------------------------------------
! NOTE: $d_{i,j}$ are already normalized, IMPORTANT to how the program !
! will compute the integrals.                                          !
! Basis sets are given with the program source code                    !
!-----------------------------------------------------------------------

CALL BASIS_MAX(basis_set, na, max_bf, max_prim, max_zeta)              ! know the maximum values to allocate the memory propertly

ALLOCATE(a_num(na), n_bf_pa(na), zeta(max_zeta,max_prim), d(max_zeta,max_prim), n_pri_bf(max_zeta))

CALL BASIS_AT(basis_set, na, atom, bf, a_num, n_bf_pa, n_pri_bf, max_bf, max_prim, max_zeta, zeta, d)

ALLOCATE(n_pri_bf_aux(bf))                                             ! Cleaning memory
                                                                       !        ⬇️
DO i = 1, bf                                                           !        ⬇️
  n_pri_bf_aux(i) = n_pri_bf(i)                                        !        ⬇️ 
ENDDO                                                                  !        ⬇️ 
                                                                       !        ⬇️ 
DEALLOCATE(n_pri_bf); ALLOCATE(n_pri_bf(bf))                           !        ⬇️ 
                                                                       !        ⬇️ 
n_pri_bf(:) = n_pri_bf_aux(:); DEALLOCATE(n_pri_bf_aux)                !  Cleaned it

!-----------------------------------------------------------------------
!********        DISTANCE BETWEEN ATOMS AND BASIS SETS          ********
!-----------------------------------------------------------------------

ALLOCATE(xyz_bf(3,bf), dis_a(na,na), dis2_a(na,na), dis_bf(bf,bf), dis2_bf(bf,bf))

CALL DISTANCES(na, bf, n_bf_pa, xyz, xyz_bf, dis_a, dis2_a, dis_bf, dis2_bf)

!-----------------------------------------------------------------------
!********            BORN-OPPENHEIMER APPROXIMATION             ********
!-----------------------------------------------------------------------

ALLOCATE(Ennmat(na,na))

CALL NNBOI(na, dis_a, a_num, Ennmat, Enn)

!-----------------------------------------------------------------------
!********                ELECTRONS AND ORBITALS                 ********
!-----------------------------------------------------------------------

n_ele = 0
DO i = 1, na                        !electrons as nuclei charge per atom
  n_ele = n_ele + a_num(i)  
ENDDO
n_ele = n_ele - q                         !do not forget the over charge

n_occ = n_ele / 2                     !how many occuped orbitals we have
n_vir = bf - n_occ                    !how many virtual orbitals we have

closed_shell = .TRUE.
IF (MODULO(n_ele,2) .NE. 0) closed_shell = .FALSE.

!-----------------------------------------------------------------------
!********                       INTEGRALS                       ********
!-----------------------------------------------------------------------

ALLOCATE(S(bf,bf), T(bf,bf), V(bf,bf), TwoEleInt(bf,bf,bf,bf))

CALL INTEGRALS(bf, na, a_num, n_pri_bf, xyz, xyz_bf, dis2_bf, max_prim, max_zeta, zeta, d, S, T, V, TwoEleInt)

!-----------------------------------------------------------------------
!********                     THEORY LEVEL                      ********
!-----------------------------------------------------------------------

SELECTCASE(method)
  CASE ('RHF')
      ALLOCATE(C(bf,bf), e(bf))
      CALL RHF(bf, n_occ, S, T, V, TwoEleInt, E_ele, C, e, step)
  CASE ('UHF')
      !CALL UHF(bf, n_occ, S, T, V, TwoEleInt, E_ele)
  CASE ('MP2')
      CALL MP2(bf, n_occ, S, T, V, TwoEleInt, E_ele, E_cmp2, step)
      E_ele = E_ele + E_cmp2
  CASE ('DFT')
      !CALL DFT(bf, n_occ, S, T, V, TwoEleInt, exchf, corrf, E_ele)
ENDSELECT
!                   ya mero cerramos el chiringuito/changarro xddxdxxdxd

!-----------------------------------------------------------------------
!********                 WRITE THE OUTPUT FILE                 ********
!-----------------------------------------------------------------------
! If all above was done correctly the work program was done,           !
!                                                          DANKE SCHÖN !

CALL WRITER(na, bf, atom, xyz, step, threshold, Enn, E_ele)

                                ENDPROGRAM nika
!                                                                 'ora sí ya cerramos el changarro/chiringuito bien /eichidi/ xdxddx
!                                                    ah nu ma' grax por leer el código si has llegado hasta acá mejor vete de fiesta
!                                                                    es más io te invito un Gin and Tonic y si no tomas pues chtm :b
!                                                                                                                   ah nocierto xdxd
!                                                                                               este programa no tiene Easter Eggs Ü
!-----------------------------------------------------------------------------------------------------------------------------------
