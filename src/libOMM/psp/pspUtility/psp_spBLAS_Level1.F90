MODULE psp_spBLAS_Level1
  use pspVariable
  use pspBasicTool
  use pspListTool

  ! This module contains sequential sparse BLAS

#ifdef MPI
  include 'mpif.h'
#endif

  private

  !**** PARAMS ************************************!

  integer, parameter :: dp=selected_real_kind(15,300)

  complex(dp), parameter :: cmplx_1=(1.0_dp,0.0_dp)
  complex(dp), parameter :: cmplx_i=(0.0_dp,1.0_dp)
  complex(dp), parameter :: cmplx_0=(0.0_dp,0.0_dp)

  !**** INTERFACES ********************************!
  interface psp_sst_DOT
     module procedure psp_sst_dDOT
     module procedure psp_sst_zDOT
  end interface psp_sst_DOT

  interface psp_sst_DOTC
     module procedure psp_sst_dDOTC
     module procedure psp_sst_zDOTC
  end interface psp_sst_DOTC

  public :: psp_sst_DOT
  public :: psp_sst_DOTC

contains

  subroutine psp_sst_dDOT(N,idxa,vala,idxb,valb,sol,numa,numb)
    ! compute the dot product of two vectors: sol=a^{T}b
    ! a is stored in a sparse vector format (idxa(1:numa),vala(1:numa))
    ! b is stored in a sparse vector format (idxb(1:numb),valb(1:numb))
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: N ! length of vectors
    integer, intent(in) :: idxa(:), idxb(:)
    real(dp), intent(in) :: vala(:), valb(:)

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: sol

    !**** OPTIONAL ********************************!

    integer, optional :: numa, numb

    !**** LOCAL ***********************************!
    integer :: ib, posa, posb, crt
    real(dp) :: va, vb

    !**********************************************!

    if (.NOT. present(numa)) then
       numa=size(vala)
    end if
    if (.NOT. present(numb)) then
       numb=size(valb)
    end if

    sol=0.0_dp
    crt=1
    do ib=1,numb
       posb=idxb(ib)
       vb=valb(ib)
       do while(crt<=numa)
          va=vala(crt)
          posa=idxa(crt)
          if (posb>posa) then
             ! next position in a
             crt=crt+1
          else if (posb==posa) then
             ! add
             sol=sol+va*vb
             crt=crt+1
             exit
          else
             exit
          end if
       end do
    end do

  end subroutine psp_sst_dDOT

  subroutine psp_sst_zDOT(N,idxa,vala,idxb,valb,sol,numa,numb)
    ! compute the dot product of two vectors: sol=a^{T}b
    ! a is stored in a sparse vector format (idxa(1:numa),vala(1:numa))
    ! b is stored in a sparse vector format (idxb(1:numb),valb(1:numb))
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: N ! length of vectors
    integer, intent(in) :: idxa(:), idxb(:)
    complex(dp), intent(in) :: vala(:), valb(:)

    !**** INOUT ***********************************!

    complex(dp), intent(inout) :: sol

    !**** OPTIONAL ********************************!

    integer, optional :: numa, numb

    !**** LOCAL ***********************************!
    integer :: ib, posa, posb, crt
    complex(dp) :: va, vb

    !**********************************************!

    if (.NOT. present(numa)) then
       numa=size(vala)
    end if
    if (.NOT. present(numb)) then
       numb=size(valb)
    end if

    sol=cmplx_0
    crt=1
    do ib=1,numb
       posb=idxb(ib)
       vb=valb(ib)
       do while(crt<=numa)
          va=vala(crt)
          posa=idxa(crt)
          if (posb>posa) then
             ! next position in a
             crt=crt+1
          else if (posb==posa) then
             ! add
             sol=sol+va*vb
             crt=crt+1
             exit
          else
             exit
          end if
       end do
    end do

  end subroutine psp_sst_zDOT

  subroutine psp_sst_dDOTC(N,idxa,vala,idxb,valb,sol,numa,numb)
    ! compute the dot product of two vectors: sol=a^{C}b
    ! a is stored in a sparse vector format (idxa(1:numa),vala(1:numa))
    ! b is stored in a sparse vector format (idxb(1:numb),valb(1:numb))    
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: N ! length of vectors
    integer, intent(in) :: idxa(:), idxb(:)
    real(dp), intent(in) :: vala(:), valb(:)

    !**** INOUT ***********************************!

    real(dp), intent(inout) :: sol

    !**** OPTIONAL ********************************!

    integer, optional :: numa, numb

    !**** LOCAL ***********************************!
    integer :: ib, posa, posb, crt
    real(dp) :: va, vb

    !**********************************************!

    if (.NOT. present(numa)) then
       numa=size(vala)
    end if
    if (.NOT. present(numb)) then
       numb=size(valb)
    end if

    sol=0.0_dp
    crt=1
    do ib=1,numb
       posb=idxb(ib)
       vb=valb(ib)
       do while(crt<=numa)
          va=vala(crt)
          posa=idxa(crt)
          if (posb>posa) then
             ! next position in a
             crt=crt+1
          else if (posb==posa) then
             ! add
             sol=sol+va*vb
             crt=crt+1
             exit
          else
             exit
          end if
       end do
    end do

  end subroutine psp_sst_dDOTC

  subroutine psp_sst_zDOTC(N,idxa,vala,idxb,valb,sol,numa,numb)
    ! compute the dot product of two vectors: sol=a^{C}b
    ! a is stored in a sparse vector format (idxa(1:numa),vala(1:numa))
    ! b is stored in a sparse vector format (idxb(1:numb),valb(1:numb))
    implicit none

    !**** INPUT ***********************************!

    integer, intent(in) :: N ! length of vectors
    integer, intent(in) :: idxa(:), idxb(:)
    complex(dp), intent(in) :: vala(:), valb(:)

    !**** INOUT ***********************************!

    complex(dp), intent(inout) :: sol

    !**** OPTIONAL ********************************!

    integer, optional :: numa, numb

    !**** LOCAL ***********************************!
    integer :: ib, posa, posb, crt
    complex(dp) :: va, vb

    !**********************************************!

    if (.NOT. present(numa)) then
       numa=size(vala)
    end if
    if (.NOT. present(numb)) then
       numb=size(valb)
    end if

    sol=cmplx_0
    crt=1
    do ib=1,numb
       posb=idxb(ib)
       vb=valb(ib)
       do while(crt<=numa)
          va=CONJG(vala(crt))
          posa=idxa(crt)
          if (posb>posa) then
             ! next position in a
             crt=crt+1
          else if (posb==posa) then
             ! add
             sol=sol+va*vb
             crt=crt+1
             exit
          else
             exit
          end if
       end do
    end do

  end subroutine psp_sst_zDOTC

END MODULE psp_spBLAS_Level1
