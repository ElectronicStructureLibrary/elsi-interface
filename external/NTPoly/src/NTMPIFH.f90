

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module wraps the MPI include statement because on certain platforms
!> just writing "USE MPI" does not work.
MODULE NTMPIModule

  IMPLICIT NONE

  INCLUDE "mpif.h"

  PUBLIC

END MODULE NTMPIModule
