MODULE pspListTool
  use pspVariable
  use pspBasicTool

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

  interface psp_spm2list
     ! create a list given a sparse matrix
     module procedure psp_dspm2list
     module procedure psp_zspm2list
  end interface psp_spm2list

  interface psp_spm2list_shift
     ! create a list given a sparse matrix
     module procedure psp_dspm2list_shift
     module procedure psp_zspm2list_shift
  end interface psp_spm2list_shift


  interface psp_spm2lists_shift
     ! create a list given a sparse matrix
     module procedure psp_dspm2lists_shift
     module procedure psp_zspm2lists_shift
  end interface psp_spm2lists_shift

  interface psp_list2spm
     ! create a list to a sparse matrix
     module procedure psp_dlist2spm
     module procedure psp_zlist2spm
  end interface psp_list2spm

  interface psp_list_create_mat
     ! create a list given a matrix
     module procedure psp_dlist_create_mat
     module procedure psp_zlist_create_mat
  end interface psp_list_create_mat

  interface psp_list_combine_listMat
     ! combine two lists, or one list one sparse matrix
     module procedure psp_dlist_combine_listMat
     module procedure psp_zlist_combine_listMat
  end interface psp_list_combine_listMat

  interface psp_list_combine_listList
     module procedure psp_dlist_combine_listList
     module procedure psp_zlist_combine_listList
  end interface psp_list_combine_listList

  interface psp_list_getLast
     ! get the last element
     module procedure psp_dlist_getLast
     module procedure psp_zlist_getLast
  end interface psp_list_getLast

  interface psp_list_print
     ! print the entire list
     module procedure psp_dlist_print
     module procedure psp_zlist_print
  end interface psp_list_print

  interface die
     module procedure die
  end interface die

  public :: psp_spm2list
  public :: psp_spm2list_shift
  public :: psp_spm2lists_shift
  public :: psp_list2spm
  public :: psp_list_create_mat
  public :: psp_list_combine_listList
  public :: psp_list_combine_listMat
  public :: psp_list_getLast
  public :: psp_list_print
  public :: psp_ddlist_combine_listMat

contains

  subroutine psp_dlist_combine_listList(m,n,alpha,listA,lenA,beta,listB,lenB,lastElem)
    ! This code essentially computes A = alpha*A+beta*B
    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: m, n
    real(dp), intent(in) :: alpha, beta

    !**** INOUT ***********************************!
    type(dList), pointer, intent(inout)     :: listA, listB, lastElem
    integer, intent(inout)    :: lenA, lenB

    !**** LOCAL ***********************************!
    integer :: iA, jA, iB, jB, cnt, crt, posA, posB
    real(dp) :: vA, vB
    type(dList), pointer  :: elem, preElem, elemB
    type(dNodeData)             :: elemData, elemAdd

    !**********************************************!
    ! multiply data
    if (alpha/=1.0_dp) then
       elem => listA
       do while( associated(elem) )
          elem%data%val=alpha*elem%data%val
          elem => list_next(elem)
       enddo
    end if

    ! begin main routine
    crt=1
    elem => listA
    preElem => listA
    elemB => listB
    do cnt=1,lenB
       iB = elemB%data%row_ind
       jB = elemB%data%col_ind
       vB = elemB%data%val
       do while(.true.)
          if ( crt<=lenA) then
             elemData = list_get_data(elem)
             iA = elemData%row_ind
             jA = elemData%col_ind
             vA = elemData%val
             posA = iA+jA*m
             posB = iB+jB*m;
             if (posB<posA) then
                elemAdd%row_ind=iB
                elemAdd%col_ind=jB
                elemAdd%val=beta*vB
                if (crt==1) then
                   call list_insert_head(listA,elemAdd)
                   preElem => listA
                else
                   call list_insert(preElem,elemAdd)
                   preElem => list_next(preElem)
                end if
                lenA=lenA+1
                crt=crt+1
                exit
             else if (posB==posA) then
                !if (beta==1.0_dp) then
                !   elem%data%val = elem%data%val + vB
                !else
                elem%data%val = elem%data%val + beta*vB
                !end if
                if (crt>1) preElem => list_next(preElem)
                elem => list_next(elem)
                elemB => list_next(elemB)
                exit
             else
                if (crt>1) preElem => list_next(preElem)
                elem => list_next(elem)
                crt = crt+1
             end if
          end if
          if ( crt>lenA ) then
             ! add the rest
             elemAdd%row_ind=iB
             elemAdd%col_ind=jB
             elemAdd%val=beta*vB
             call list_insert(preElem,elemAdd)
             lenA=lenA+1
             crt=crt+1
             preElem => list_next(preElem)
             exit
          end if
       end do
       elemB => list_next(elemB)
    end do
    lastElem => preElem

  end subroutine psp_dlist_combine_listList

  subroutine psp_zlist_combine_listList(m,n,alpha,listA,lenA,beta,listB,lenB,lastElem)
    ! This code essentially computes A = alpha*A+beta*B
    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: m, n
    complex(dp), intent(in) :: alpha, beta

    !**** INOUT ***********************************!
    type(zList), pointer, intent(inout)     :: listA, listB, lastElem
    integer, intent(inout)    :: lenA, lenB

    !**** LOCAL ***********************************!
    integer :: iA, jA, iB, jB, cnt, crt, posA, posB
    complex(dp) :: vA, vB
    type(zList), pointer  :: elem, preElem, elemB
    type(zNodeData)             :: elemData, elemAdd

    !**********************************************!
    ! multiply data
    if (alpha/=cmplx_1) then
       elem => listA
       do while( associated(elem) )
          elem%data%val=alpha*elem%data%val
          elem => list_next(elem)
       enddo
    end if

    ! begin main routine
    crt=1
    elem => listA
    preElem => listA
    elemB => listB
    do cnt=1,lenB
       iB = elemB%data%row_ind
       jB = elemB%data%col_ind
       vB = elemB%data%val
       do while(.true.)
          if ( crt<=lenA) then
             elemData = list_get_data(elem)
             iA = elemData%row_ind
             jA = elemData%col_ind
             vA = elemData%val
             posA = iA+jA*m
             posB = iB+jB*m;
             if (posB<posA) then
                elemAdd%row_ind=iB
                elemAdd%col_ind=jB
                elemAdd%val=beta*vB
                if (crt==1) then
                   call list_insert_head(listA,elemAdd)
                   preElem => listA
                else
                   call list_insert(preElem,elemAdd)
                   preElem => list_next(preElem)
                end if
                lenA=lenA+1
                crt=crt+1
                exit
             else if (posB==posA) then
                !if (beta==1.0_dp) then
                !   elem%data%val = elem%data%val + vB
                !else
                elem%data%val = elem%data%val + beta*vB
                !end if
                if (crt>1) preElem => list_next(preElem)
                elem => list_next(elem)
                elemB => list_next(elemB)
                exit
             else
                if (crt>1) preElem => list_next(preElem)
                elem => list_next(elem)
                crt = crt+1
             end if
          end if
          if ( crt>lenA ) then
             ! add the rest
             elemAdd%row_ind=iB
             elemAdd%col_ind=jB
             elemAdd%val=beta*vB
             call list_insert(preElem,elemAdd)
             lenA=lenA+1
             crt=crt+1
             preElem => list_next(preElem)
             exit
          end if
       end do
       elemB => list_next(elemB)
    end do
    lastElem => preElem

  end subroutine psp_zlist_combine_listList

  subroutine psp_dlist_create_mat(listA,lenA,beta,idxi,idxj,val,fmt,numAddInB,lastElem)
    ! This code create a list with a given sparse matrix (idxi,idxj,beta*val).
    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: idxi(:), idxj(:)
    real(dp), intent(in) :: val(:), beta
    character(3), intent(in) :: fmt ! format of B
    integer, intent(in) :: numAddInB
    ! the first numAddInB elements in the sparse matrix are added in list

    !**** INOUT ***********************************!
    type(dList), pointer, intent(inout)   :: listA, lastElem
    integer, intent(inout)    :: lenA

    !**** LOCAL ***********************************!
    integer :: iB, jB, cnt, numAdded
    real(dp) :: vB
    type(dList), pointer  :: elem
    type(dNodeData)             :: elemAdd

    !**********************************************!
    if (numAddInB>0) then
       ! initialization
       lenA=0
       !if ( associated(listA) ) listA=>null()!call list_destroy(listA)!segmentation fault if use list_destroy
       ! begin main routine
       numAdded=0
       if (fmt=='csc') then
          do jB=1,size(idxj)-1
             do cnt=idxj(jB),idxj(jB+1)-1
                iB = idxi(cnt)
                vB = val(cnt)
                elemAdd%row_ind=iB
                elemAdd%col_ind=jB
                elemAdd%val=beta*vB
                if (numAdded==0) then
                   call list_create(listA,elemAdd)
                   elem => listA
                else
                   call list_insert(elem,elemAdd)
                   elem => list_next(elem)
                end if
                lenA=lenA+1
                numAdded=numAdded+1
                if (numAdded>=numAddInB) exit
             end do
             if (numAdded>=numAddInB) exit
          end do
       else
          do cnt=1,size(val)
             iB = idxi(cnt)
             jB = idxj(cnt)
             vB = val(cnt)
             elemAdd%row_ind=iB
             elemAdd%col_ind=jB
             elemAdd%val=beta*vB
             if (numAdded==0) then
                call list_create(listA,elemAdd)
                elem => listA
             else
                call list_insert(elem,elemAdd)
                elem => list_next(elem)
             end if
             lenA=lenA+1
             numAdded=numAdded+1
             if (numAdded>=numAddInB) exit
          end do
       end if
       lastElem => elem
    else
       lastElem=>listA
    end if

  end subroutine psp_dlist_create_mat

  subroutine psp_zlist_create_mat(listA,lenA,beta,idxi,idxj,val,fmt,numAddInB,lastElem)
    ! This code create a list with a given sparse matrix (idxi,idxj,beta*val).
    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: idxi(:), idxj(:)
    complex(dp), intent(in) :: val(:), beta
    character(3), intent(in) :: fmt ! format of B
    integer, intent(in) :: numAddInB
    ! the first numAddInB elements in the sparse matrix are added in list

    !**** INOUT ***********************************!
    type(zList), pointer, intent(inout)   :: listA, lastElem
    integer, intent(inout)    :: lenA

    !**** LOCAL ***********************************!
    integer :: iB, jB, cnt, numAdded
    complex(dp) :: vB
    type(zList), pointer  :: elem
    type(zNodeData)             :: elemAdd

    !**********************************************!
    if (numAddInB>0) then
       ! initialization
       lenA=0
       !if ( associated(listA) ) listA=>null()!call list_destroy(listA)!segmentation fault if use list_destroy
       ! begin main routine
       numAdded=0
       if (fmt=='csc') then
          do jB=1,size(idxj)-1
             do cnt=idxj(jB),idxj(jB+1)-1
                iB = idxi(cnt)
                vB = val(cnt)
                elemAdd%row_ind=iB
                elemAdd%col_ind=jB
                elemAdd%val=beta*vB
                if (numAdded==0) then
                   call list_create(listA,elemAdd)
                   elem => listA
                else
                   call list_insert(elem,elemAdd)
                   elem => list_next(elem)
                end if
                lenA=lenA+1
                numAdded=numAdded+1
                if (numAdded>=numAddInB) exit
             end do
             if (numAdded>=numAddInB) exit
          end do
       else
          do cnt=1,size(val)
             iB = idxi(cnt)
             jB = idxj(cnt)
             vB = val(cnt)
             elemAdd%row_ind=iB
             elemAdd%col_ind=jB
             elemAdd%val=beta*vB
             if (numAdded==0) then
                call list_create(listA,elemAdd)
                elem => listA
             else
                call list_insert(elem,elemAdd)
                elem => list_next(elem)
             end if
             lenA=lenA+1
             numAdded=numAdded+1
             if (numAdded>=numAddInB) exit
          end do
       end if
       lastElem => elem
    else
       lastElem=>listA
    end if

  end subroutine psp_zlist_create_mat

  subroutine psp_dlist_combine_listMat(m,n,alpha,listA,lenA,beta,idxi,idxj,val,fmt,numAddInB,lastElem)
    ! This code essentially computes A = alpha*A+beta*B,
    ! if A is stored in a list format.
    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: idxi(:), idxj(:), m, n
    real(dp), intent(in) :: val(:), alpha, beta
    character(3), intent(in) :: fmt ! format of B
    integer, intent(in) :: numAddInB
    ! the first numAddInB elements in the sparse matrix are added in list

    !**** INOUT ***********************************!
    type(dList), pointer, intent(inout)   :: listA, lastElem
    integer, intent(inout)    :: lenA

    !**** LOCAL ***********************************!
    integer :: iA, jA, iB, jB, cnt, crt, posA, posB, numAdded
    real(dp) :: vA, vB
    type(dList), pointer  :: elem, preElem
    type(dNodeData)             :: elemData, elemAdd

    !**********************************************!

    ! multiply data
    if (alpha/=1.0_dp) then
       numAdded=0
       elem => listA
       do while( associated(elem) .and. numAdded<numAddInB )
          numAdded=numAdded+1
          elem%data%val=alpha*elem%data%val
          elem => list_next(elem)
       enddo
    end if

    ! begin main routine
    crt=1
    numAdded=0
    elem => listA
    preElem => listA
    if (numAddInB>0) then
       if (fmt=='csc') then
          do jB=1,n
             do cnt=idxj(jB),idxj(jB+1)-1
                iB = idxi(cnt)
                vB = val(cnt)
                numAdded=numAdded+1
                do while(.true.)
                   if ( crt<=lenA) then
                      elemData = list_get_data(elem)
                      iA = elemData%row_ind
                      jA = elemData%col_ind
                      vA = elemData%val
                      posA = iA+jA*m
                      posB = iB+jB*m;
                      if (posB<posA) then
                         elemAdd%row_ind=iB
                         elemAdd%col_ind=jB
                         elemAdd%val=beta*vB
                         if (crt==1) then
                            call list_insert_head(listA,elemAdd)
                            preElem => listA
                         else
                            call list_insert(preElem,elemAdd)
                            preElem => list_next(preElem)
                         end if
                         lenA=lenA+1
                         crt=crt+1
                         exit
                      else if (posB==posA) then
                         elem%data%val = elem%data%val + beta*vB
                         if (crt>1) preElem => list_next(preElem)
                         elem => list_next(elem)
                         crt=crt+1
                         exit
                      else
                         if (crt>1) preElem => list_next(preElem)
                         elem => list_next(elem)
                         crt = crt+1
                      end if
                   end if
                   if ( crt>lenA ) then
                      ! add the rest
                      elemAdd%row_ind=iB
                      elemAdd%col_ind=jB
                      elemAdd%val=beta*vB
                      call list_insert(preElem,elemAdd)
                      lenA=lenA+1
                      crt=crt+1
                      preElem => list_next(preElem)
                      exit
                   end if
                end do
                if (numAdded>=numAddInB) exit
             end do
             if (numAdded>=numAddInB) exit
          end do
       else
          do cnt=1,size(val)
             iB = idxi(cnt)
             jB = idxj(cnt)
             vB = val(cnt)
             numAdded=numAdded+1
             do while(.true.)
                if ( crt<=lenA) then
                   elemData = list_get_data(elem)
                   iA = elemData%row_ind
                   jA = elemData%col_ind
                   vA = elemData%val
                   posA = iA+jA*m
                   posB = iB+jB*m;
                   if (posB<posA) then
                      elemAdd%row_ind=iB
                      elemAdd%col_ind=jB
                      elemAdd%val=beta*vB
                      if (crt==1) then
                         call list_insert_head(listA,elemAdd)
                         preElem => listA
                      else
                         call list_insert(preElem,elemAdd)
                         preElem => list_next(preElem)
                      end if
                      lenA=lenA+1
                      crt=crt+1
                      exit
                   else if (posB==posA) then
                      elem%data%val = elem%data%val + beta*vB
                      if (crt>1) preElem => list_next(preElem)
                      elem => list_next(elem)
                      crt=crt+1
                      exit
                   else
                      if (crt>1) preElem => list_next(preElem)
                      elem => list_next(elem)
                      crt = crt+1
                   end if
                end if
                if ( crt>lenA ) then
                   ! add the rest
                   elemAdd%row_ind=iB
                   elemAdd%col_ind=jB
                   elemAdd%val=beta*vB
                   call list_insert(preElem,elemAdd)
                   lenA=lenA+1
                   crt=crt+1
                   preElem => list_next(preElem)
                   exit
                end if
             end do
             if (numAdded>=numAddInB) exit
          end do
       end if
    end if
    if (crt>lenA) then
       lastElem => preElem
    else
       lastElem => elem
    end if
    call psp_list_getLast(listA,lastElem)

  end subroutine psp_dlist_combine_listMat

  subroutine psp_zlist_combine_listMat(m,n,alpha,listA,lenA,beta,idxi,idxj,val,fmt,numAddInB,lastElem)
    ! This code essentially computes A = alpha*A+beta*B
    implicit none

    !**** INPUT ***********************************!
    integer, intent(in) :: idxi(:), idxj(:), m, n
    complex(dp), intent(in) :: val(:), alpha, beta
    character(3), intent(in) :: fmt ! format of B
    integer, intent(in) :: numAddInB
    ! the first numAddInB elements in the sparse matrix are added in list

    !**** INOUT ***********************************!
    type(zList), pointer, intent(inout)   :: listA, lastElem
    integer, intent(inout)    :: lenA

    !**** LOCAL ***********************************!
    integer :: iA, jA, iB, jB, cnt, crt, posA, posB, numAdded
    complex(dp) :: vA, vB
    type(zList), pointer  :: elem, preElem
    type(zNodeData)             :: elemData, elemAdd

    !**********************************************!

    ! multiply data
    if (alpha/=cmplx_1) then
       elem => listA
       do while( associated(elem) )
          elem%data%val=alpha*elem%data%val
          elem => list_next(elem)
       enddo
    end if

    ! begin main routine
    crt=1
    numAdded=0
    elem => listA
    preElem => listA
    if ( numAddInB>0) then
       if (fmt=='csc') then
          do jB=1,n
             do cnt=idxj(jB),idxj(jB+1)-1
                iB = idxi(cnt)
                vB = val(cnt)
                numAdded=numAdded+1
                do while(.true.)
                   if ( crt<=lenA) then
                      elemData = list_get_data(elem)
                      iA = elemData%row_ind
                      jA = elemData%col_ind
                      vA = elemData%val
                      posA = iA+jA*m
                      posB = iB+jB*m
                      if (posB<posA) then
                         elemAdd%row_ind=iB
                         elemAdd%col_ind=jB
                         elemAdd%val=beta*vB
                         if (crt==1) then
                            call list_insert_head(listA,elemAdd)
                            preElem => listA
                         else
                            call list_insert(preElem,elemAdd)
                            preElem => list_next(preElem)
                         end if
                         lenA=lenA+1
                         crt=crt+1
                         !numAdded=numAdded+1
                         exit
                      else if (posB==posA) then
                         elem%data%val = elem%data%val + beta*vB
                         if (crt>1) preElem => list_next(preElem)
                         elem => list_next(elem)
                         crt=crt+1
                         !numAdded=numAdded+1
                         exit
                      else
                         if (crt>1) preElem => list_next(preElem)
                         elem => list_next(elem)
                         crt = crt+1
                      end if
                   end if
                   if ( crt>lenA ) then
                      ! add the rest
                      elemAdd%row_ind=iB
                      elemAdd%col_ind=jB
                      elemAdd%val=beta*vB
                      call list_insert(preElem,elemAdd)
                      lenA=lenA+1
                      crt=crt+1
                      preElem => list_next(preElem)
                      !numAdded=numAdded+1
                      exit
                   end if
                end do
                if (numAdded>=numAddInB) exit
             end do
             if (numAdded>=numAddInB) exit
          end do
       else
          do cnt=1,size(val)
             iB = idxi(cnt)
             jB = idxj(cnt)
             vB = val(cnt)
             numAdded=numAdded+1
             do while(.true.)
                if ( crt<=lenA) then
                   elemData = list_get_data(elem)
                   iA = elemData%row_ind
                   jA = elemData%col_ind
                   vA = elemData%val
                   posA = iA+jA*m
                   posB = iB+jB*m;
                   if (posB<posA) then
                      elemAdd%row_ind=iB
                      elemAdd%col_ind=jB
                      elemAdd%val=beta*vB
                      if (crt==1) then
                         call list_insert_head(listA,elemAdd)
                         preElem => listA
                      else
                         call list_insert(preElem,elemAdd)
                         preElem => list_next(preElem)
                      end if
                      lenA=lenA+1
                      crt=crt+1
                      !numAdded=numAdded+1
                      exit
                   else if (posB==posA) then
                      elem%data%val = elem%data%val + beta*vB
                      if (crt>1) preElem => list_next(preElem)
                      elem => list_next(elem)
                      crt=crt+1
                      !numAdded=numAdded+1
                      exit
                   else
                      if (crt>1) preElem => list_next(preElem)
                      elem => list_next(elem)
                      crt = crt+1
                   end if
                end if
                if ( crt>lenA ) then
                   ! add the rest
                   elemAdd%row_ind=iB
                   elemAdd%col_ind=jB
                   elemAdd%val=beta*vB
                   call list_insert(preElem,elemAdd)
                   lenA=lenA+1
                   crt=crt+1
                   preElem => list_next(preElem)
                   !numAdded=numAdded+1
                   exit
                end if
             end do
             if (numAdded>=numAddInB) exit
          end do
       end if
    end if
    if (crt>lenA) then
       lastElem => preElem
    else
       lastElem => elem
    end if
    call psp_list_getLast(listA,lastElem)


  end subroutine psp_zlist_combine_listMat

  subroutine psp_dspm2list_shift(m,n,idx1,idx2,val,nnz,fmt,list,last,disp)
    ! This code convert a sparse matrix into a list format.
    ! If disp/=[0,], then the matrix is treated as a submatrix and
    ! its position is shifted disp(1) down and disp(2) to the right.

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt! storage format of sparse matrices, 'coo' or 'csc'
    real(dp), intent(in) :: val(:)
    integer, intent(in) :: m, n
    integer, intent(in) :: idx1(:), idx2(:)

    !**** INOUT ***********************************!
    type(dList), pointer, intent(inout)  :: list, last
    integer, intent(inout) :: nnz

    !**** OPTIONAL ***********************************!
    integer, dimension(2), optional :: disp

    !**** INTERNAL ********************************!
    type(dNodeData)      :: node
    integer :: cnt, cnt2, i, j
    real(dp) :: v
    type(dList), pointer :: elem

    !**********************************************!

    if (.NOT. present(disp)) then
       disp(1)=0
       disp(2)=0
    end if

    nnz = size(val)
    if (nnz>0) then
       select case (fmt)
       case ('coo')
          node%val = val(1)
          node%col_ind = idx2(1)+disp(2)
          node%row_ind = idx1(1)+disp(1)
          call list_create( list, node )
          elem => list
          if (nnz>1) then
             do cnt=2,nnz
                node%val = val(cnt)
                node%col_ind = idx2(cnt)+disp(2)
                node%row_ind = idx1(cnt)+disp(1)
                call list_insert( elem, node )
                elem=>list_next(elem)
             end do
          end if
          last=>elem
       case ('csc')
          do j=1,n
             do cnt=idx2(j),idx2(j+1)-1
                i=idx1(cnt)
                v=val(cnt)
                ! add (i,j,v) to list
                if (cnt==1) then
                   node%val = v
                   node%col_ind = j+disp(2)
                   node%row_ind = i+disp(1)
                   call list_create( list, node )
                   elem => list
                else
                   node%val = v
                   node%col_ind = j+disp(2)
                   node%row_ind = i+disp(1)
                   call list_insert( elem, node )
                   elem => list_next(elem)
                end if
             end do
          end do
          last=>elem
       end select
    end if

  end subroutine psp_dspm2list_shift

  subroutine psp_zspm2list_shift(m,n,idx1,idx2,val,nnz,fmt,list,last,disp)
    ! This code convert a sparse matrix into a list format.
    ! If disp/=[0,], then the matrix is treated as a submatrix and
    ! its position is shifted disp(1) down and disp(2) to the right.

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt! storage format of sparse matrices, 'coo' or 'csc'
    complex(dp), intent(in) :: val(:)
    integer, intent(in) :: m, n
    integer, intent(in) :: idx1(:), idx2(:)

    !**** INOUT ***********************************!
    type(zList), pointer, intent(inout)  :: list, last
    integer, intent(inout) :: nnz

    !**** OPTIONAL ***********************************!
    integer, dimension(2), optional :: disp

    !**** INTERNAL ********************************!
    type(zNodeData)      :: node
    integer :: cnt, cnt2, i, j
    complex(dp) :: v
    type(zList), pointer :: elem

    !**********************************************!

    if (.NOT. present(disp)) then
       disp(1)=0
       disp(2)=0
    end if

    nnz = size(val)
    if (nnz>0) then
       select case (fmt)
       case ('coo')
          node%val = val(1)
          node%col_ind = idx2(1)+disp(2)
          node%row_ind = idx1(1)+disp(1)
          call list_create( list, node )
          elem => list
          if (nnz>1) then
             do cnt=2,nnz
                node%val = val(cnt)
                node%col_ind = idx2(cnt)+disp(2)
                node%row_ind = idx1(cnt)+disp(1)
                call list_insert( elem, node )
                elem=>list_next(elem)
             end do
          end if
          last=>elem
       case ('csc')
          do j=1,n
             do cnt=idx2(j),idx2(j+1)-1
                i=idx1(cnt)
                v=val(cnt)
                ! add (i,j,v) to list
                if (cnt==1) then
                   node%val = v
                   node%col_ind = j+disp(2)
                   node%row_ind = i+disp(1)
                   call list_create( list, node )
                   elem => list
                else
                   node%val = v
                   node%col_ind = j+disp(2)
                   node%row_ind = i+disp(1)
                   call list_insert( elem, node )
                   elem => list_next(elem)
                end if
             end do
          end do
          last=>elem
       end select
    end if

  end subroutine psp_zspm2list_shift

  subroutine psp_dspm2lists_shift(m,n,idx1,idx2,val,nnz,fmt,listArraySt,listArrayEd,listCreated,disp)
    ! This code convert a sparse matrix into a list format.
    ! If disp/=[0,], then the matrix is treated as a submatrix and
    ! its position is shifted disp(1) down and disp(2) to the right.

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt! storage format of sparse matrices, 'coo' or 'csc'
    real(dp), intent(in) :: val(:)
    integer, intent(in) :: m, n
    integer, intent(in) :: idx1(:), idx2(:)

    !**** INOUT ***********************************!
    integer, intent(inout) :: nnz
    type(dListPtrArray), dimension(n), intent(inout) :: listArraySt, listArrayEd
    logical, dimension(n) :: listCreated

    !**** OPTIONAL ***********************************!
    integer, dimension(2), optional :: disp

    !**** INTERNAL ********************************!
    type(dNodeData)      :: node
    integer :: cnt, cnt2, i, j
    real(dp) :: v

    !**********************************************!

    if (.NOT. present(disp)) then
       disp(1)=0
       disp(2)=0
    end if

    do cnt=1,n
       listCreated(cnt)=.false.
    end do

    nnz = size(val)
    if (nnz>0) then
       select case (fmt)
       case ('coo')
          do cnt=1,nnz
             node%val = val(cnt)
             node%col_ind = idx2(cnt)+disp(2)
             node%row_ind = idx1(cnt)+disp(1)
             if (.not.listCreated(idx2(cnt))) then
                listCreated(idx2(cnt))=.true.
                call list_create( listArraySt(idx2(cnt))%ptr, node )
                listArrayEd(idx2(cnt))%ptr => listArraySt(idx2(cnt))%ptr
             else
                call list_insert(listArrayEd(idx2(cnt))%ptr, node )
                listArrayEd(idx2(cnt))%ptr => list_next(listArrayEd(idx2(cnt))%ptr)
             end if
          end do

       case ('csc')
          do j=1,n
             do cnt=idx2(j),idx2(j+1)-1
                i=idx1(cnt)
                v=val(cnt)
                ! add (i,j,v) to list
                if (.not.listCreated(j)) then
                   listCreated(j)=.true.
                   node%val = v
                   node%col_ind = j+disp(2)
                   node%row_ind = i+disp(1)
                   call list_create( listArraySt(j)%ptr, node )
                   listArrayEd(j)%ptr => listArraySt(j)%ptr
                else
                   node%val = v
                   node%col_ind = j+disp(2)
                   node%row_ind = i+disp(1)
                   call list_insert(listArrayEd(j)%ptr, node )
                   listArrayEd(j)%ptr => list_next(listArrayEd(j)%ptr)
                end if
             end do
          end do
       end select
    end if

  end subroutine psp_dspm2lists_shift

  subroutine psp_zspm2lists_shift(m,n,idx1,idx2,val,nnz,fmt,listArraySt,listArrayEd,listCreated,disp)
    ! This code convert a sparse matrix into a list format.
    ! If disp/=[0,], then the matrix is treated as a submatrix and
    ! its position is shifted disp(1) down and disp(2) to the right.

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt! storage format of sparse matrices, 'coo' or 'csc'
    complex(dp), intent(in) :: val(:)
    integer, intent(in) :: m, n
    integer, intent(in) :: idx1(:), idx2(:)

    !**** INOUT ***********************************!
    integer, intent(inout) :: nnz
    type(zListPtrArray), dimension(n), intent(inout) :: listArraySt, listArrayEd
    logical, dimension(n) :: listCreated

    !**** OPTIONAL ***********************************!
    integer, dimension(2), optional :: disp

    !**** INTERNAL ********************************!
    type(zNodeData)      :: node
    integer :: cnt, cnt2, i, j
    complex(dp) :: v

    !**********************************************!

    if (.NOT. present(disp)) then
       disp(1)=0
       disp(2)=0
    end if

    do cnt=1,n
       listCreated(cnt)=.false.
    end do

    nnz = size(val)
    if (nnz>0) then
       select case (fmt)
       case ('coo')
          do cnt=1,nnz
             node%val = val(cnt)
             node%col_ind = idx2(cnt)+disp(2)
             node%row_ind = idx1(cnt)+disp(1)
             if (.not.listCreated(idx2(cnt))) then
                listCreated(idx2(cnt))=.true.
                call list_create( listArraySt(idx2(cnt))%ptr, node )
                listArrayEd(idx2(cnt))%ptr => listArraySt(idx2(cnt))%ptr
             else
                call list_insert(listArrayEd(idx2(cnt))%ptr, node )
                listArrayEd(idx2(cnt))%ptr => list_next(listArrayEd(idx2(cnt))%ptr)
             end if
          end do

       case ('csc')
          do j=1,n
             do cnt=idx2(j),idx2(j+1)-1
                i=idx1(cnt)
                v=val(cnt)
                ! add (i,j,v) to list
                if (.not.listCreated(j)) then
                   listCreated(j)=.true.
                   node%val = v
                   node%col_ind = j+disp(2)
                   node%row_ind = i+disp(1)
                   call list_create( listArraySt(j)%ptr, node )
                   listArrayEd(j)%ptr => listArraySt(j)%ptr
                else
                   node%val = v
                   node%col_ind = j+disp(2)
                   node%row_ind = i+disp(1)
                   call list_insert(listArrayEd(j)%ptr, node )
                   listArrayEd(j)%ptr => list_next(listArrayEd(j)%ptr)
                end if
             end do
          end do
       end select
    end if

  end subroutine psp_zspm2lists_shift


  subroutine psp_dspm2list(m,n,idx1,idx2,val,nnz,fmt,list)

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt! storage format of sparse matrices, 'coo' or 'csc'
    real(dp), intent(in) :: val(:)
    integer, intent(in) :: m, n
    integer, intent(in) :: idx1(:), idx2(:)

    !**** INOUT ***********************************!
    type(dList), pointer, intent(inout)  :: list
    integer, intent(inout) :: nnz

    !**** INTERNAL ********************************!
    type(dNodeData)      :: node
    integer :: cnt, cnt2, i, j
    real(dp) :: v
    type(dList), pointer :: elem

    !**********************************************!

    nnz = size(val)
    if (nnz>0) then
       select case (fmt)
       case ('coo')
          node%val = val(1)
          node%col_ind = idx2(1)
          node%row_ind = idx1(1)
          call list_create( list, node )
          elem => list
          if (nnz>1) then
             do cnt=2,nnz
                node%val = val(cnt)
                node%col_ind = idx2(cnt)
                node%row_ind = idx1(cnt)
                call list_insert( elem, node )
                elem=>list_next(elem)
             end do
          end if
       case ('csc')
          do j=1,n
             do cnt=idx2(j),idx2(j+1)-1
                i=idx1(cnt)
                v=val(cnt)
                ! add (i,j,v) to list
                if (cnt==1) then
                   node%val = v
                   node%col_ind = j
                   node%row_ind = i
                   call list_create( list, node )
                   elem => list
                else
                   node%val = v
                   node%col_ind = j
                   node%row_ind = i
                   call list_insert( elem, node )
                   elem => list_next(elem)
                end if
             end do
          end do
       end select
    end if

  end subroutine psp_dspm2list

  subroutine psp_zspm2list(m,n,idx1,idx2,val,nnz,fmt,list)

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt! storage format of sparse matrices, 'coo' or 'csc'
    complex(dp), intent(in) :: val(:)
    integer, intent(in) :: m, n
    integer, intent(in) :: idx1(:), idx2(:)

    !**** INOUT ***********************************!
    type(zList), pointer, intent(inout)  :: list
    integer, intent(inout) :: nnz

    !**** INTERNAL ********************************!
    type(zNodeData)      :: node
    integer :: cnt, cnt2, i, j
    complex(dp) :: v
    type(zList), pointer :: elem

    !**********************************************!

    nnz = size(val)
    if (nnz>0) then
       select case (fmt)
       case ('coo')
          node%val = val(1)
          node%col_ind = idx2(1)
          node%row_ind = idx1(1)
          call list_create( list, node )
          elem => list
          if (nnz>1) then
             do cnt=2,nnz
                node%val = val(cnt)
                node%col_ind = idx2(cnt)
                node%row_ind = idx1(cnt)
                call list_insert( elem, node )
                elem=>list_next(elem)
             end do
          end if
       case ('csc')
          do j=1,n
             do cnt=idx2(j),idx2(j+1)-1
                i=idx1(cnt)
                v=val(cnt)
                ! add (i,j,v) to list
                if (cnt==1) then
                   node%val = v
                   node%col_ind = j
                   node%row_ind = i
                   call list_create( list, node )
                   elem => list
                else
                   node%val = v
                   node%col_ind = j
                   node%row_ind = i
                   call list_insert( elem, node )
                   elem => list_next(elem)
                end if
             end do
          end do
       end select
    end if

  end subroutine psp_zspm2list


  subroutine psp_dlist2spm(m,n,idx1,idx2,val,fmt,list,nnz,ifCount)

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt! storage format of sparse matrices, 'coo' or 'csc'
    integer, intent(in) :: m, n
    type(dList), pointer, intent(in)  :: list
    logical, intent(in)    :: ifCount ! .true. count the length of list
    ! .false. use the nnz provided

    !**** INOUT ***********************************!
    real(dp), allocatable, intent(inout) :: val(:)
    integer, allocatable, intent(inout) :: idx1(:), idx2(:)
    integer, intent(inout) :: nnz

    !**** INTERNAL ********************************!
    type(dList), pointer  :: elem
    type(dNodeData)             :: data
    integer :: cnt

    !**********************************************!

    if (ifCount.eqv..true.) nnz = list_count(list)
    if (allocated(val)) deallocate(val)
    if (allocated(idx1)) deallocate(idx1)
    if (allocated(idx2)) deallocate(idx2)
    allocate(val(nnz))
    allocate(idx1(nnz))
    allocate(idx2(nnz))

    if (nnz>0) then
       cnt=0
       elem => list
       !do while( associated(elem) )
       do while( cnt<nnz )
          cnt=cnt+1
          data = list_get_data(elem)
          idx1(cnt)=data%row_ind
          idx2(cnt)=data%col_ind
          val(cnt)=data%val
          elem => list_next(elem)
       end do
    end if
    if (fmt=='csc') call psp_sst_coo2csc(m,n,nnz,idx2)


  end subroutine psp_dlist2spm

  subroutine psp_zlist2spm(m,n,idx1,idx2,val,fmt,list,nnz,ifCount)

    !**** INPUT ***********************************!
    character(3), intent(in) ::   fmt! storage format of sparse matrices, 'coo' or 'csc'
    integer, intent(in) :: m, n
    type(zList), pointer, intent(in)  :: list
    logical, intent(in)    :: ifCount ! .true. count the length of list
    ! .false. use the nnz provided

    !**** INOUT ***********************************!
    complex(dp), allocatable, intent(inout) :: val(:)
    integer, allocatable, intent(inout) :: idx1(:), idx2(:)
    integer, intent(inout) :: nnz

    !**** INTERNAL ********************************!
    type(zList), pointer  :: elem
    type(zNodeData)             :: data
    integer :: cnt

    !**********************************************!

    if (ifCount.eqv..true.) nnz = list_count(list)
    if (allocated(val)) deallocate(val)
    if (allocated(idx1)) deallocate(idx1)
    if (allocated(idx2)) deallocate(idx2)
    allocate(val(nnz))
    allocate(idx1(nnz))
    allocate(idx2(nnz))

    if (nnz>0) then
       cnt=0
       elem => list
       !do while( associated(elem) )
       do while( cnt<nnz )
          cnt=cnt+1
          data = list_get_data(elem)
          idx1(cnt)=data%row_ind
          idx2(cnt)=data%col_ind
          val(cnt)=data%val
          elem => list_next(elem)
       end do
    end if
    if (fmt=='csc') call psp_sst_coo2csc(m,n,nnz,idx2)

  end subroutine psp_zlist2spm

  subroutine psp_dlist_getLast(list,lastElem)
    type(dList), pointer, intent(in)  :: list

    type(dList), pointer, intent(inout)  :: lastElem

    type(dList), pointer  :: preElem

    lastElem => list
    do while( associated(lastElem) )
       preElem => lastElem
       lastElem => list_next(lastElem)
    end do
    lastElem => preElem

  end subroutine psp_dlist_getLast

  subroutine psp_zlist_getLast(list,lastElem)
    type(zList), pointer, intent(in)  :: list

    type(zList), pointer, intent(inout)  :: lastElem

    type(zList), pointer  :: preElem

    lastElem => list
    do while( associated(lastElem) )
       preElem => lastElem
       lastElem => list_next(lastElem)
    end do
    lastElem => preElem

  end subroutine psp_zlist_getLast


  subroutine psp_dlist_print( text, list, len )
    character(len=*)                 :: text
    integer, optional :: len
    type(dList), pointer  :: list

    type(dList), pointer  :: elem
    type(dNodeData)             :: data
    integer                          :: i

    write(*,*) text
    if (present(len)) print *, 'length is ', len
    !
    ! Loop over the list
    !
    i = 0
    elem => list
    do while( associated(elem) )
       data = list_get_data(elem)
       write(*,*) i, data%row_ind, ',   ', data%col_ind, ',   ', data%val
       elem => list_next(elem)
       i = i + 1
    enddo

  end subroutine psp_dlist_print

  subroutine psp_zlist_print( text, list, len )
    character(len=*)                 :: text
    integer, optional :: len
    type(zList), pointer  :: list

    type(zList), pointer  :: elem
    type(zNodeData)             :: data
    integer                          :: i

    write(*,*) text
    if (present(len)) print *, 'length is ', len
    !
    ! Loop over the list
    !
    i = 0
    elem => list
    do while( associated(elem) )
       data = list_get_data(elem)
       write(*,*) i, data%row_ind, ',   ', data%col_ind, ',   ', data%val
       elem => list_next(elem)
       i = i + 1
    enddo

  end subroutine psp_zlist_print

  subroutine die(message)
    implicit none

    !**** INPUT ***********************************!

    character(*), intent(in), optional :: message

    !**** INTERNAL ********************************!

    logical, save :: log_start=.false.

    integer :: log_unit

    !**********************************************!

    if (log_start) then
       open(newunit=log_unit,file='MatrixSwitch.log',position='append')
    else
       open(newunit=log_unit,file='MatrixSwitch.log',status='replace')
       log_start=.true.
    end if
    write(log_unit,'(a)'), 'FATAL ERROR in matrix_switch!'
    write(log_unit,'(a)'), message
    close(log_unit)
    stop

  end subroutine die


END MODULE pspListTool
