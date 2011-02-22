PROGRAM ACCQUAREL
  USE setup_tools ; USE case_parameters ; USE data_parameters ; USE basis_parameters ; USE scf_parameters
  USE random
  IMPLICIT NONE
  DOUBLE PRECISION :: STARTTIME,ENDTIME

  CALL CPU_TIME(STARTTIME)

  CALL INIT_RANDOM_SEED
  
  WRITE(*,*)' *** ACCQUAREL using A.S.P.I.C. *** '
  !$ WRITE(*,*)' (parallel version compiled by an OpenMP-compliant implementation)'
! Determination of the name of the setup file
  CALL SETUP_FILENAME
! Case setup
  CALL SETUP_CASE
! Approximation setup
  CALL SETUP_APPROXIMATION
! Formalism/model setup
  CALL SETUP_FORMALISM
! SCF algorithms setup
  CALL SETUP_SCF
! Molecular system setup
  CALL SETUP_DATA
! Discretization basis setup
  CALL SETUP_BASIS
! Computation
  IF (RELATIVISTIC) THEN
     CALL DRIVER_relativistic
  ELSE
     IF (APPROXIMATION==1) THEN
        CALL DRIVER_boson_star
     ELSE
        CALL DRIVER_nonrelativistic
     END IF
  END IF

  CALL CPU_TIME(ENDTIME)
  WRITE(*,*)' ' ; WRITE(*,*)' Elapsed CPU time is ',ENDTIME-STARTTIME,' seconds.'
END PROGRAM ACCQUAREL
