module pspListType

  use psp_dList, only: dList   => LINKED_LIST,   &
       dNodeData          => LIST_DATA,     &
       dlist_create   => list_create,   &
       dlist_insert   => list_insert,   &
       dlist_insert_head   => list_insert_head,   &
       dlist_count    => list_count,    &
       dlist_next     => list_next,     &
       dlist_get_data => list_get_data, &
       dlist_put_data => list_put_data, &
       dlist_delete_element   => list_delete_element, &
       dlist_destroy  => list_destroy

  use psp_zList, only: zList   => LINKED_LIST,   &
       zNodeData            => LIST_DATA,     &
       zlist_create     => list_create,   &
       zlist_insert     => list_insert,   &
       zlist_insert_head   => list_insert_head,   &
       zlist_count      => list_count,    &
       zlist_next       => list_next,     &
       zlist_get_data   => list_get_data, &
       zlist_put_data   => list_put_data, &
       zlist_delete_element     => list_delete_element, &
       zlist_destroy    => list_destroy
  private

  type dListPtrArray
     type(dList), pointer :: ptr
  end type dListPtrArray

  type zListPtrArray
     type(zList), pointer :: ptr
  end type zListPtrArray

  public :: list_create, list_destroy, list_get_data, list_put_data, &
       list_insert, list_insert_head, list_next, list_delete_element, list_count
  public :: dList, zList
  public :: dNodeData, zNodeData
  public :: dListPtrArray, zListPtrArray

  interface list_create
     module procedure dlist_create
     module procedure zlist_create
  end interface list_create

  interface list_count
     module procedure dlist_count
     module procedure zlist_count
  end interface list_count

  interface list_destroy
     module procedure dlist_destroy
     module procedure zlist_destroy
  end interface list_destroy

  interface list_insert
     module procedure dlist_insert
     module procedure zlist_insert
  end interface list_insert

  interface list_insert_head
     module procedure dlist_insert_head
     module procedure zlist_insert_head
  end interface list_insert_head

  interface list_next
     module procedure dlist_next
     module procedure zlist_next
  end interface list_next

  interface list_get_data
     module procedure dlist_get_data
     module procedure zlist_get_data
  end interface list_get_data

  interface list_put_data
     module procedure dlist_put_data
     module procedure zlist_put_data
  end interface list_put_data

  interface list_delete_element
     module procedure dlist_delete_element
     module procedure zlist_delete_element
  end interface list_delete_element

end module pspListType
