module energy_grid

  use ace_header,       only: RrrData
  use constants,        only: MAX_LINE_LEN
  use global
  use list_header,      only: ListReal
  use output,           only: write_message

contains

!===============================================================================
! UNIONIZED_GRID creates a single unionized energy grid combined from each
! nuclide of each material. Right now, the grid for each nuclide is added into a
! linked list one at a time with an effective insertion sort. Could be done with
! a hash for all energy points and then a quicksort at the end (what hash
! function to use?)
!===============================================================================

  subroutine unionized_grid()

    integer :: i ! index in nuclides array
    type(ListReal), pointer :: list => null()
    type(Nuclide),  pointer :: nuc => null()

    message = "Creating unionized energy grid..."
    call write_message(5)

    ! Add grid points for each nuclide in the problem
    do i = 1, n_nuclides_total
      nuc => nuclides(i)
      call add_grid_points(list, nuc % energy)
    end do

    ! Set size of unionized energy grid 
    n_grid = list % size() 

    ! create allocated array from linked list
    allocate(e_grid(n_grid))
    do i = 1, n_grid
      e_grid(i) = list % get_item(i)
    end do

    ! delete linked list and dictionary
    call list % clear()
    deallocate(list)

    ! Set pointers to unionized energy grid for each nuclide
    call grid_pointers()

    write (*,*) "Size of unionized grid is", n_grid ! delete

  end subroutine unionized_grid

!===============================================================================
! ADD_GRID_POINTS adds energy points from the 'energy' array into a linked list
! of points already stored from previous arrays.
!===============================================================================

  subroutine add_grid_points(list, energy)

    type(ListReal), pointer :: list
    real(8), intent(in) :: energy(:)

    integer :: i       ! index in energy array
    integer :: n       ! size of energy array
    integer :: current ! current index 
    real(8) :: E       ! actual energy value

    i = 1
    n = size(energy)

    ! If the original list is empty, we need to allocate the first element and
    ! store first energy point
    if (.not. associated(list)) then
      allocate(list)
      do i = 1, n
        call list % append(energy(i))
      end do
      return
    end if

    ! Set current index to beginning of the list 
    current = 1

    do while (i <= n)
      E = energy(i)

      ! If we've reached the end of the grid energy list, add the remaining
      ! energy points to the end
      if (current > list % size()) then
        ! Finish remaining energies
        do while (i <= n)
          call list % append(energy(i))
          i = i + 1
        end do
        exit
      end if

      if (E < list % get_item(current)) then

        ! Insert new energy in this position
        call list % insert(current, E)

        ! Advance index in linked list and in new energy grid
        i = i + 1
        current = current + 1

      elseif (E == list % get_item(current)) then
        ! Found the exact same energy, no need to store duplicates so just
        ! skip and move to next index
        i = i + 1
        current = current + 1
      else
        current = current + 1
      end if

    end do

  end subroutine add_grid_points

!===============================================================================
! GRID_POINTERS creates an array of pointers (ints) for each nuclide to link
! each point on the nuclide energy grid to one on the unionized energy grid
! Updated to work with clustered data
! Index of the left nuclide gridpoint given if union gridpoint is between two.
!===============================================================================

  subroutine grid_pointers()

    integer :: i            ! loop index for nuclides
    integer :: j            ! loop index for union energy grid
    integer :: index_e      ! index of the nuclide energy grid (right index)
    real(8) :: union_energy ! energy on union grid
    real(8) :: energy       ! energy on nuclide grid
    integer :: offset_c     ! offset into the codebook (clustering only)
    integer :: offset_fast  ! offset into the fast region (clustering only)
    type(Nuclide), pointer :: nuc => null()
    type(RrrData), pointer :: rrr => null()

    do i = 1, n_nuclides_total
      nuc => nuclides(i)
      allocate(nuc % grid_index(n_grid))

      index_e = 1
      energy = nuc % energy(index_e)

      if (nuc % rrr_cluster) then
        ! Within the clustered region the codebook is used for the grid index
        ! Discontiguities of energy clusters cause more points to be in
        ! the nuclide energy grid (in the RRR) than there are clusters, which
        ! means offsets must be used in the indexing.
        rrr => nuc % rrr_data
        offset_c = rrr % i_low
        offset_fast = rrr % n_clust - (rrr % i_high - rrr % i_low)
        !write (*,*) nuc % name, rrr % i_low, rrr % i_high, rrr % n_clust, &
        !  nuc % n_grid
        do j = 1, n_grid
          union_energy = e_grid(j)
          if (union_energy >= energy .and. index_e < nuc % n_grid) then
            index_e = index_e + 1
            energy = nuc % energy(index_e)
          end if
          if (index_e <= rrr % i_low) then
            ! thermal region, no offset
            nuc % grid_index(j) = index_e - 1
          else if (index_e > rrr % i_high) then
            ! fast region, offset
            nuc % grid_index(j) = index_e + offset_fast - 1
          else
            ! RRR region, index from codebook plus offset
            nuc % grid_index(j) = rrr % codebook(index_e - offset_c) + &
              offset_c - 1
          end if
        end do
      else
        ! No clustering. Use normal indexing.
        do j = 1, n_grid
          union_energy = e_grid(j)
          if (union_energy >= energy .and. index_e < nuc % n_grid) then
            index_e = index_e + 1
            energy = nuc % energy(index_e)
          end if
          nuc % grid_index(j) = index_e - 1
        end do
      end if
      !write (*,*) '---------------------------------'
      !write (*,*) nuc % name
      !write (*,'(8e14.6)') nuc % energy(1:200)
      !write (*,*) '---------------------------------'
      !write (*,'(8e14.6)') e_grid(1:200)
      !write (*,'(8i14)') nuc % grid_index(1:200)
      !write (*,*) '---------------------------------'
    end do



  end subroutine grid_pointers

end module energy_grid
