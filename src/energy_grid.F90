module energy_grid

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
!===============================================================================

  subroutine grid_pointers()

    integer :: i            ! loop index for nuclides
    integer :: j            ! loop index for nuclide energy grid
    integer :: index_e      ! index on union energy grid
    real(8) :: union_energy ! energy on union grid
    real(8) :: energy       ! energy on nuclide grid
    type(Nuclide), pointer :: nuc => null()

    do i = 1, n_nuclides_total
      nuc => nuclides(i)
      allocate(nuc % grid_index(n_grid))

      index_e = 1
      energy = nuc % energy(index_e)

      do j = 1, n_grid
        union_energy = e_grid(j)
        if (union_energy >= energy .and. index_e < nuc % n_grid) then
          index_e = index_e + 1
          energy = nuc % energy(index_e)
        end if
        nuc % grid_index(j) = index_e - 1
      end do
    end do

  end subroutine grid_pointers

!===============================================================================
! CASCADING_GRID creates augmented nuclide energy grids so that XS lookups can
! be done using fractional cascading
!===============================================================================

  subroutine cascading_grid()

    integer :: i            ! index in nuclides array
    integer :: j            ! index in augmented energy grid
    integer :: i_nuc        ! current index in nuclide energy grid
    integer :: i_aug        ! current index in next augmented energy grid
    real(8) :: E_nuc        ! energy in nuclide grid
    real(8) :: E_aug        ! energy in next augmented grid
    type(Nuclide), pointer :: nuc => null() ! pointer to nuclide i
    type(Nuclide), pointer :: nuc_next => null() ! pointer to nuclide i + 1

    message = "Creating cascading energy grid..."
    call write_message(5)

    do i = n_nuclides_total, 1, -1
      nuc => nuclides(i)

      ! The last augmented list is just the original list
      if (i == n_nuclides_total) then
        nuc % n_aug_grid = nuc % n_grid
        allocate(nuc % aug_energy(nuc % n_aug_grid))
        allocate(nuc % aug_pointers(2, nuc % n_aug_grid))
        nuc % aug_energy = nuc % energy

        ! fill the array with pointer pairs
        nuc % aug_pointers(1, :) = (/(j, j = 1, nuc % n_aug_grid)/)
        nuc % aug_pointers(2, :) = 1

      ! Create the augmented nuclide energy grid and pointer arrays
      else
        nuc_next => nuclides(i + 1)
        nuc % n_aug_grid = nuc % n_grid + nuc_next % n_aug_grid / 2
        allocate(nuc % aug_energy(nuc % n_aug_grid))
        allocate(nuc % aug_pointers(2, nuc % n_aug_grid))

        ! Set indices to the beginning of the nuclide energy arrays
        j = 1
        i_nuc = 1
        ! Set current index of next augmented array to 2 since we will only
        ! add every other element to the current augmented array
        i_aug = 2

        do while (j <= nuc % n_aug_grid)

          ! If we've reached the end of nuclide energy grid, add the
          ! remaining energy points from the next augmented grid to the end
          if (i_nuc > nuc % n_grid) then
            do while (j <= nuc % n_aug_grid)
              nuc % aug_energy(j) = nuc_next % aug_energy(i_aug)
              ! Set nuclide grid pointer out of bounds since energy is beyond it
              nuc % aug_pointers(1, j) = i_nuc
              nuc % aug_pointers(2, j) = i_aug
              i_aug = i_aug + 2
              j = j + 1
            end do
            exit
          end if

          ! If we've reached the end of the next augmented  grid, add the
          ! remaining energy points from the nuclide grid to the end   
          if (i_aug > nuc_next % n_aug_grid) then
            do while (j <= nuc % n_aug_grid)
              nuc % aug_energy(j) = nuc % energy(i_nuc)
              ! Set next augmented grid pointer to out of bounds
              nuc % aug_pointers(1, j) = i_nuc
              nuc % aug_pointers(2, j) = i_aug
              i_nuc = i_nuc + 1
              j = j + 1
            end do
            exit
          end if

          E_nuc = nuc % energy(i_nuc)
          E_aug = nuc_next % aug_energy(i_aug)

          ! Set the nuclide grid and next augmented grid pointers
          nuc % aug_pointers(1, j) = i_nuc
          nuc % aug_pointers(2, j) = i_aug

          ! Add the smaller of the nuclide and next augmented grid energies
          ! and advance that index
          if(E_nuc <= E_aug) then
            nuc % aug_energy(j) = E_nuc
            i_nuc = i_nuc + 1
            j = j + 1
          else
            nuc % aug_energy(j) = E_aug
            i_aug = i_aug + 2
            j = j + 1
          end if

        end do

      end if

    end do

  end subroutine cascading_grid

end module energy_grid
