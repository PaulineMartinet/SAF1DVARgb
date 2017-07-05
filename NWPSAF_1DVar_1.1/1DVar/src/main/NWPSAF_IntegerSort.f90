SUBROUTINE NWPSAF_IntegerSort( &
  key,   & ! in
  index)   ! inout

!
!    This software was developed within the context of 
!    the EUMETSAT Satellite Application Facility on 
!    Numerical Weather Prediction (NWP SAF), under the 
!    Cooperation Agreement dated 25 November 1998, between 
!    EUMETSAT and the Met Office, UK, by one or more partners 
!    within the NWP SAF. The partners in the NWP SAF are 
!    the Met Office, ECMWF, KNMI and MeteoFrance.
! 
!    Copyright 2004, EUMETSAT, All Rights Reserved.
!
! Description:
!   Generates a index array pointing to the elements of the array 'key'
!   in increasing order
!
! Method:
!   The heap sort invented by J.W.J.Williams is used.
!   A description of the method can be found in 'Numerical Recipes'
!   The group information array is used to allow easy sorting on several
!   parameters of different types.
!
! Inputs: 
!   key : An array of character strings, to be sorted
!
! Input/Output: 
!   index : An integer array pointing to the sorted items.
!
! Owner: Manager of operational observation processing and assimilation
!
! History:
!
! Version Date     Comment
! ------- ----     -------
! 1.0     10/06/02 First version, based on Ops_IntegerSortQuick. A. Collard.
!
! Code Description:
!   Language:       Fortran 90 
!   Software Standards: GTDP 8
!
! End of header -------------------------------------------------------

! Used modules:


IMPLICIT NONE

! Arguments:
INTEGER, INTENT(IN) :: key(:)
INTEGER, POINTER :: index(:)

! Local variables:
INTEGER :: n     ! The number of items
INTEGER :: head  ! heaps are tree structures: head and child refer
INTEGER :: child ! to related items within the tree 
INTEGER :: j          
INTEGER :: dum   ! used to swap index items

!Could put in an optional mask
n=SIZE(key)
IF (.NOT.ASSOCIATED(Index)) ALLOCATE(Index(n))
DO j=1,n
  Index(j)=j
END DO

! Do heapsort: Create the heap...
makeheap : DO j=n/2,1,-1
  head=j
  sift1 : DO
     ! find the largest out of the head and its two children...
    child=head*2
    IF (child>n) EXIT sift1
    IF (child<n) THEN
      IF (key(Index(child+1))>key(Index(child))) child=child+1
    END IF
     ! if the head is the largest, then sift is done...
    IF (key(Index(head))>=key(Index(child))) EXIT sift1
     ! otherwise swap to put the largest child at the head,
     ! and prepare to repeat the procedure for the head in its new
     ! subordinate position.
    dum=Index(child)
    Index(child)=Index(head)
    Index(head)=dum
    head=child
  END DO sift1
END DO makeheap
 ! Retire heads of the heap, which are the largest, and
 ! stack them at the end of the array.
retire : DO j=n,2,-1
  dum=Index(1)
  Index(1)=Index(j)
  Index(j)=dum
  head=1
   ! second sift is similar to first...
  sift2: DO
    child=head*2
    IF (child>(j-1)) EXIT sift2
    IF (child<(j-1)) THEN
      IF (key(Index(child+1))>key(Index(child))) child=child+1
    END IF
    IF (key(Index(head))>=key(Index(child))) EXIT sift2
    dum=Index(child)
    Index(child)=Index(head)
    Index(head)=dum
    head=child
  END DO sift2  
END DO retire

END SUBROUTINE NWPSAF_IntegerSort
