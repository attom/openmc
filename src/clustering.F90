module clustering

  use ace_header,        only: Nuclide, RrrData, Reaction
  use clustering_kmeans, only: perform_kms, kms_uniform_clust_cen
  use constants
  use error,             only: fatal_error
  use global
  use search,            only: binary_search
  use string,            only: to_str

  implicit none

contains

!===============================================================================
! CLUSTER_ONE_NUCLIDE clusters the cross sections of one nuclide in the resolved
! resonance region
!===============================================================================

  subroutine cluster_one_nuclide(i_nuclide)

    integer, intent(in) :: i_nuclide       ! index of nuclide to be clustered
    real(8), allocatable :: observations(:,:) ! cross sections to be clustered
    real(8), allocatable :: clust_cen(:,:) ! centroids of clusters
    integer :: n_rrr                       ! size of RRR grid before clustering
    integer :: n_rxn                       ! number of reactions to cluster
    logical err_flag ! false if clustering was successful
    real(8) :: SMALL = 1E-11_8 ! a small value so log calls are robust
    type(Nuclide), pointer :: nuc => null() ! pointer to nuclide
    type(RrrData), pointer :: rrr => null() ! pointer to nuclide's clust data
    type(Reaction), pointer :: rxn => null() ! pointer to nuclide's rxn data
    integer :: i,j
    integer :: threshold ! threshold index for reactions

    nuc => nuclides(i_nuclide)
    ! Do not cluster light elements for now
    if (.not. clustering_on .or. nuc % zaid <= 11000) return
    !
    write (*,*) "Clustering nuclide ", nuc % name
    !write (*,*) "Using ", n_clust_glob, &
    !  " clusters and ", n_group_glob, " groups."

    ! Allocate RRR data for nuclide
    nuc % rrr_cluster = .true.
    allocate(nuc % rrr_data)
    rrr => nuc % rrr_data

    ! Determine range of RRR (in MeV) from URR, if available
    rrr % e_low = 4 / 1E6_8 ! Make a parameter later
    rrr % e_high = 20 ! Make a parameter later
    if (nuc % urr_present) then
      rrr % e_high = nuc % urr_data % energy(1)
    end if
    rrr % i_low = find_grid_index(nuc % energy, &
      rrr % e_low)
    rrr % i_high = find_grid_index(nuc % energy, &
      rrr % e_high)
    if (nuc % urr_present) then
      if (nuc % energy(rrr % i_high) >= nuc % urr_data % energy(1)) then
        rrr % i_high = rrr % i_high - 1
      end if
    end if

    ! Do not cluster in the fast energy region, as defined by rxn thresholds.
    do i = 1, nuc % n_reaction
      rxn => nuc % reactions(i)
      threshold = rxn % threshold
      if (threshold /= 1 .and. rxn % MT < N_GAMMA) then
        rrr % i_high = min(rrr % i_high, threshold - 1)
      end if
    end do
    if (rrr % i_high <= rrr % i_low) then 
      nuc % rrr_cluster = .false.
      return
    end if
    rrr % e_low = nuc % energy(rrr % i_low)
    rrr % e_high = nuc % energy(rrr % i_high)
    !
    !write (*,*) "RRR goes from ", rrr % e_low, " MeV to ", &
    !  rrr % e_high, " MeV."
    write (*,*) "On the nuclide grid, this is ", &
      nuc % energy(rrr % i_low), " MeV to ", &
      nuc % energy(rrr % i_high), " MeV."
    !write (*,*) "(indices ", rrr % i_low, " - ", &
    !  rrr % i_high, ")"

    ! Copy cross sections in RRR to observations buffer
    ! We want to preserve cross sections on a log-log scale
    ! We inclue total, elastic, and maybe fission
    n_rxn = 2
    if (nuc % fissionable) n_rxn = 3
    n_rrr = rrr % i_high - rrr % i_low + 1
    allocate(observations(n_rxn, n_rrr))
    observations(1, :) = log10(nuc % total(rrr%i_low : rrr%i_high) + SMALL)
    observations(2, :) = log10(nuc % elastic(rrr%i_low : rrr%i_high) + SMALL)
    if (nuc % fissionable) then
      observations(3, :) = log10(nuc % fission(rrr%i_low : rrr%i_high) + SMALL)
    end if
    !
    !write (*,*) minval(observations), maxval(observations)

    ! Perform k-means clustering
    rrr % n_clust = n_clust_glob
    call kms_uniform_clust_cen(observations, rrr % n_clust, clust_cen)
    call perform_kms(observations, clust_max_it, clust_tol, clust_cen, &
      rrr % codebook, rrr % L2_err, rrr % Linf_err, err_flag)
    if (err_flag) then
      message = 'K-means failed to produce a valid clustering'
      call fatal_error()
    end if
    !
    !do i = 1,size(clust_cen, 2)
    !  write(*, '(100e9.2)') (10**clust_cen(j, i), j=1,n_rxn)
    !end do
    write (*,*) "Error is", rrr % L2_err, "/", rrr % Linf_err

    call write_clustering(i_nuclide, observations, clust_cen)

    ! Overwrite pointwise cross section data with clustered data
    call apply_clustering_to_all_xs(i_nuclide)
    !nuc % n_grid = n_xs

    !do i = 1, nuc % n_reaction
    !  write (*,*) 'MT', nuc % reactions(i) % MT, &
    !  'with size', size(nuc % reactions(i) % sigma), &
    !  'and threshold', nuc % reactions(i) % threshold, &
    !  '/', nuc % n_grid
    !end do

    ! Update indices that point to the nuclide's energy grid
    ! Reactions -> threshold (integer)
    ! Energy grid (update and thin).
    ! Energy grid fast offset
    write (*,*) '----------------------------------'

  end subroutine cluster_one_nuclide

  subroutine write_clustering(i_nuclide, observations, clust_cen)

    integer, intent(in) :: i_nuclide       ! index of nuclide to be clustered
    real(8), allocatable :: observations(:,:) ! cross sections to be clustered
    real(8), allocatable :: clust_cen(:,:) ! centroids of clusters
    real(8) :: SMALL = 1E-11_8 ! a small value so log calls are robust
    type(Nuclide), pointer :: nuc => null() ! pointer to nuclide

    nuc => nuclides(i_nuclide)
    open(11, file='obs_' // trim(adjustl(nuc % name)) // '.txt')
    open(12, file='cen_' // trim(adjustl(nuc % name)) // '.txt')
    if (nuc % fissionable) then
      write(11, '(3e12.4)') 10**observations - SMALL
      write(12, '(3e12.4)') 10**clust_cen - SMALL
    else
      write(11, '(2f12.3)') 10**observations - SMALL
      write(12, '(2f12.3)') 10**clust_cen - SMALL
    end if
    close(11)
    close(12)

  end subroutine write_clustering


!===============================================================================
! APPLY_CLUSTERING_TO_ALL_XS uses the codebook to cluster all the cross sections
! for one nuclide
!===============================================================================

  subroutine apply_clustering_to_all_xs(i_nuclide)

    integer, intent(in) :: i_nuclide ! index of nuclide to be clustered
    integer :: n_rrr   ! size of RRR grid before clustering
    integer :: n_therm ! size of xs below the RRR
    integer :: n_fast  ! size of xs above RRR (includes URR and fast)
    integer :: n_clust ! size of RRR after clustering
    type(Nuclide), pointer :: nuc => null() ! pointer to nuclide
    type(RrrData), pointer :: rrr => null() ! pointer to nuclide's clust data
    !integer :: threshold ! threshold index for reactions
    integer :: n_xs_old, n_xs_new  ! delete

    ! Determine sizes and set pointers
    nuc => nuclides(i_nuclide)
    rrr => nuc % rrr_data
    n_therm = rrr % i_low - 1
    n_fast = nuc % n_grid - rrr % i_high
    n_clust = rrr % n_clust
    n_rrr = rrr % i_high - rrr % i_low + 1

    n_xs_old = n_therm + n_rrr + n_fast ! delete
    n_xs_new = n_therm + n_clust + n_fast ! delete
    write(*,*) 'old size', n_xs_old, 'new size', n_xs_new, &
      'clusters', rrr % n_clust

    ! Apply clusteirng for total, elastic, etc. reactions for the nuclide
    call condense_one_xs(nuc % total, rrr % codebook, n_therm, n_fast, &
      n_clust, n_rrr)

  end subroutine apply_clustering_to_all_xs

!===============================================================================
! CONDENSE_ONE_XS uses the codebook to condense the RRR of one cross section.
! The energy range is split into three pieces: thermal, RRR, and fast (which
! includes URR). The thermal and fast ranges are simply copied. The RRR is
! condensed into n_clust numbers using the codebook.
!===============================================================================

  subroutine condense_one_xs(xs, codebook, n_therm, n_fast, n_clust, n_rrr)

    real(8), allocatable, intent(inout) :: xs(:) ! xs to be condensed
    integer, allocatable, intent(in) :: codebook(:) ! how to do the condensation
    integer, intent(in) :: n_rrr   ! size of RRR grid before clustering
    integer, intent(in) :: n_therm ! size of xs below the RRR
    integer, intent(in) :: n_fast  ! size of xs above RRR (includes URR and fast)
    integer, intent(in) :: n_clust ! size of RRR after clustering
    real(8), allocatable :: xs_scratch(:)  ! scratch space for one cross section
    integer, allocatable :: xs_norm(:)     ! denominator for weighted average
    integer :: wgt ! weight for the weighted cross section averaging
    integer :: n_xs_old, n_xs_new  ! xs sizes
    integer :: offset_rrr, offset_fast_old, offset_fast_new ! xs array offsets
    integer :: i, strt, endd, i_clust ! indices for arrays and loops
    real(8) :: SMALL = 1E-11_8 ! a small value so log calls are robust

    ! Determine offsets and do allocations
    offset_rrr = n_therm
    offset_fast_old = offset_rrr + n_rrr
    offset_fast_new = offset_rrr + n_clust
    n_xs_old = n_therm + n_rrr + n_fast
    n_xs_new = n_therm + n_clust + n_fast
    allocate(xs_scratch(n_xs_new))
    allocate(xs_norm(n_clust))
    xs_scratch = 0
    xs_norm = 0

    ! Copy thermal
    xs_scratch(1:n_therm) = xs(1:n_therm)

    ! Cluster RRR
    wgt = 1
    do i = 1, n_rrr
      i_clust = codebook(i)
      xs_norm(i_clust) = xs_norm(i_clust) + wgt
      xs_scratch(i_clust + offset_rrr) = xs_scratch(i_clust + offset_rrr) + &
        log(xs(i + offset_rrr)+SMALL) * wgt
    end do
    strt = offset_rrr + 1
    endd = offset_rrr + n_clust
    xs_norm = max(xs_norm, 1)
    xs_scratch(strt:endd) = exp(xs_scratch(strt:endd) / xs_norm)

    ! Copy URR and fast
    xs_scratch(offset_fast_new:n_xs_new) = xs(offset_fast_old:n_xs_old)

    ! Move (requires Fortran 2003 compiler but prevents another copy)
    deallocate(xs)
    call move_alloc(xs_scratch, xs)

    !write (*,*) 'Cross section'
    !do i = offset_rrr+1, offset_rrr + n_clust
    !  write(*,*) xs(i)
    !end do

  end subroutine condense_one_xs

!===============================================================================
! FIND_GRID_INDEX determines the index of a sorted grid closest to an input
! value
!===============================================================================

  function find_grid_index(grid, val) result (loc)

    real(8), intent(in), allocatable :: grid(:) ! sorted grid
    real(8), intent(in) :: val                  ! value to look up
    integer :: loc                              ! grid(loc) is close to val
    integer :: n_grid_local                     ! size of grid

    n_grid_local = size(grid)

    ! if particle's energy is outside of energy grid range, set to first or last
    ! index. Otherwise, do a binary search through the union energy grid.
    if (val < grid(1)) then
      loc = 1
    elseif (val > grid(n_grid_local)) then
      loc = n_grid_local
    else
      loc = binary_search(grid, n_grid_local, val)
    end if

  end function find_grid_index

end module clustering
