module clustering

  use ace_header,        only: Nuclide, RrrData, Reaction
  use clustering_kmeans, only: perform_kms, kms_uniform_clust_cen
  use constants
  use error,             only: fatal_error
  use global
  use output,            only: write_message
  use search,            only: binary_search
  use string,            only: to_str

  implicit none

contains

!===============================================================================
! CLUSTER_ALL_NUCLIDES clusters the XS in the resolved resonance region (RRR)
! for all the nuclides.
!===============================================================================

  subroutine cluster_all_nuclides()
    integer :: i_nuclide

    ! Assumes the nuclides go from 1 to n_nuclides_total.
!$omp parallel do schedule(dynamic)
    do i_nuclide = 1, n_nuclides_total
      call cluster_one_nuclide(i_nuclide)
    end do
!$omp end parallel do

  end subroutine

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
    type(Nuclide), pointer, save :: nuc => null() ! pointer to nuclide
    type(RrrData), pointer, save :: rrr => null() ! ptr to nuclide's clust data
    type(Reaction), pointer, save :: rxn => null() ! ptr to nuclide's rxn data
    integer :: threshold ! integer where the reaction begins in the energy grid
    integer :: i
!$omp threadprivate(nuc, rrr, rxn)

    nuc => nuclides(i_nuclide)
    ! Do not cluster light elements for now
    if (nuc % zaid <= 11000) return

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

    !write (*,*) "RRR goes from ", rrr % e_low, " MeV to ", &
    !  rrr % e_high, " MeV."
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

    ! Determine the number of clusters
    rrr % n_clust = n_clust_glob
    if (rrr % n_clust >= n_rrr) then
        nuc % rrr_cluster = .false.
        return
    end if

    ! Perform k-means clustering
    call kms_uniform_clust_cen(observations, rrr % n_clust, clust_cen)
    call perform_kms(observations, clust_max_it, clust_tol, clust_cen, &
      rrr % codebook, rrr % L2_err, rrr % Linf_err, err_flag)
    if (err_flag) then
      message = 'K-means failed to produce a valid clustering'
      call fatal_error()
    end if

    !do i = 1,size(clust_cen, 2)
    !  write(*, '(100e9.2)') (10**clust_cen(j, i), j=1,n_rxn)
    !end do

    ! Write to output files the observations and clusterings
    !call write_clustering(i_nuclide, observations, clust_cen)

    !do i = 1, nuc % n_reaction
    !  write (*,*) 'MT', nuc % reactions(i) % MT, &
    !  'with size', size(nuc % reactions(i) % sigma), &
    !  'and threshold', nuc % reactions(i) % threshold, &
    !  '/', nuc % n_grid
    !end do

    ! Overwrite pointwise cross section data with clustered data
    call apply_clustering_to_all_xs(i_nuclide)

    !do i = 1, nuc % n_reaction
    !  write (*,*) 'MT', nuc % reactions(i) % MT, &
    !  'with size', size(nuc % reactions(i) % sigma), &
    !  'and threshold', nuc % reactions(i) % threshold, &
    !  '/', nuc % n_grid
    !end do

    ! Write brief description
    write (*,*) "Clustering nuclide ", nuc % name
    write (*,*) "On the nuclide grid, this is ", &
      nuc % energy(rrr % i_low), " MeV to ", &
      nuc % energy(rrr % i_high), " MeV."
    write (*,*) "Error is", rrr % L2_err, "/", rrr % Linf_err
    write (*,*) '----------------------------------'

  end subroutine cluster_one_nuclide

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
    integer i
    type(Nuclide), pointer, save :: nuc => null() ! pointer to nuclide
    type(RrrData), pointer, save :: rrr => null() ! ptr to nuclide's clust data
    type(Reaction), pointer, save :: rxn => null() ! ptr to nuclide reaction dat
    integer :: n_xs_old, n_xs_new  ! delete
!$omp threadprivate(nuc, rrr, rxn)

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

    ! Re-order clusters for cache-friendliness.
    call reorder_clusters(rrr % codebook, n_therm, n_rrr, n_clust)

    ! Do the identity mapping for verification (delete)
    !n_clust = n_rrr
    !do i = 1, size(rrr % codebook)
    !  rrr%codebook(i) = i
    !end do

    ! Apply clusteirng for total, elastic, etc. xs arrays for the nuclide
    call condense_one_xs(nuc % total, rrr % codebook, n_therm, n_fast, &
      n_clust, n_rrr)
    call condense_one_xs(nuc % elastic, rrr % codebook, n_therm, n_fast, &
      n_clust, n_rrr)
    call condense_one_xs(nuc % absorption, rrr % codebook, n_therm, n_fast, &
      n_clust, n_rrr)
    call condense_one_xs(nuc % heating, rrr % codebook, n_therm, n_fast, &
      n_clust, n_rrr)
    call condense_one_xs(nuc % fission, rrr % codebook, n_therm, n_fast, &
      n_clust, n_rrr)
    call condense_one_xs(nuc % nu_fission, rrr % codebook, n_therm, n_fast, &
      n_clust, n_rrr)

    ! Apply clustering to each reaction's sigma array. The sizes will be
    ! reaction-dependent because only cross sections after the threshold are
    ! stored.
    do i = 1, nuc % n_reaction
      rxn => nuc % reactions(i)
      if (allocated(rxn % sigma) .and. rxn % threshold < rrr % i_low) then
        n_therm = rrr % i_low - rxn % threshold
        n_rrr = rrr % i_high - rrr % i_low + 1
        n_fast = nuc % n_grid - rrr % i_high
        !write (*,*) rxn % MT
        call condense_one_xs(rxn % sigma, rrr % codebook, n_therm, n_fast, &
          n_clust, n_rrr)
      else if (rxn % MT < N_GAMMA .and. &
          rxn % threshold <= rrr % i_high .and. &
          rxn % threshold >= rrr % i_low) then
        !write (*,*) rxn % MT, rxn % threshold, rrr % i_low, rrr % i_high
        message = "Assumed no thresholds in the RRR for clustering. Aborting."
        call fatal_error()
      end if
    end do

    ! Thin the nuclide energy grid. Update codebook and thresholds for thinned
    ! array.
    call thin_energy_grid(i_nuclide)

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

    ! If cross section is not allocated, skip
    if (.not. allocated(xs)) return

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

    ! Copy thermal portion, which may not exist for thresholded reactions
    if (n_therm > 0) xs_scratch(1:n_therm) = xs(1:n_therm)

    ! Cluster RRR
    ! We do a weighted *linear* averaging to preserve (sum of partials) equals
    ! total in each cluster
    wgt = 1
    do i = 1, n_rrr
      i_clust = codebook(i)
      xs_norm(i_clust) = xs_norm(i_clust) + wgt
      xs_scratch(i_clust + offset_rrr) = xs_scratch(i_clust + offset_rrr) + &
        xs(i + offset_rrr) * wgt
    end do
    strt = offset_rrr + 1
    endd = offset_rrr + n_clust
    xs_norm = max(xs_norm, 1)
    xs_scratch(strt:endd) = xs_scratch(strt:endd) / xs_norm

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
! THIN_ENERGY_GRID thins the energy grid to remove gridpoints that share
! the same codebook value
!===============================================================================
  subroutine thin_energy_grid(i_nuclide)

    integer, intent(in) :: i_nuclide ! index of nuclide to be clustered
    integer :: n_rrr   ! size of RRR grid before clustering
    integer :: n_therm ! size of xs below the RRR
    integer :: n_fast  ! size of xs above RRR (includes URR and fast)
    integer :: n_clust ! number of clusters
    integer :: n_grid_old, n_grid_new ! sizes of nuclide's energy grid
    real(8), allocatable :: e_scratch(:)  ! scratch space for energy grid
    integer, allocatable :: c_scratch(:) ! scratch for thinned grid codebook
    type(Nuclide), pointer, save :: nuc => null() ! pointer to nuclide
    type(RrrData), pointer, save :: rrr => null() ! ptr to nuclide's clust data
    type(Reaction), pointer, save :: rxn => null() ! ptr to nuclide reaction dat
    integer :: i_cur ! current location in e_scratch / c_scratch
    integer :: i_clust, i_clust_old ! indices for cluster
    integer :: offset_rrr, offset_fast ! xs array offsets
    integer :: i, strt, endd ! indices for arrays and loops
!$omp threadprivate(nuc, rrr, rxn)

    ! Determine sizes and set pointers and offsets
    nuc => nuclides(i_nuclide)
    rrr => nuc % rrr_data
    n_therm = rrr % i_low - 1
    n_fast = nuc % n_grid - rrr % i_high
    n_rrr = rrr % i_high - rrr % i_low + 1
    n_clust = rrr % n_clust
    offset_rrr = n_therm
    offset_fast = n_therm + n_rrr
    n_grid_old = size(nuc % energy)

    ! We do not upfront know how much thinning will take place,
    ! but the thinned array will be at most the size of energy.
    allocate(e_scratch(n_grid_old))
    e_scratch = ZERO

    ! Similarly, the thinned codebook will be at most the size of the old
    ! codebook.
    allocate(c_scratch(n_rrr))
    c_scratch = 0

    ! Copy thermal portion of energy range
    e_scratch(1:n_therm) = nuc % energy(1:n_therm)

    ! Thin the RRR portion of the energy range. Only need to keep energy
    ! gridpoints where the cluster number changes.
    i_clust_old = 0
    i_cur = 0
    do i = 1, n_rrr
      i_clust = rrr % codebook(i)
      if (i_clust /= i_clust_old) then
        i_clust_old = i_clust
        i_cur = i_cur + 1
        e_scratch(offset_rrr + i_cur) = nuc % energy(offset_rrr + i)
        c_scratch(i_cur) = i_clust
      end if
    end do

    ! Copy fast portion of the energy range
    strt = offset_rrr + i_cur + 1
    endd = offset_rrr + i_cur + n_fast
    e_scratch(strt:endd) =  nuc % energy(offset_fast+1 : offset_fast+n_fast)

    ! Determine new length for energy grid
    n_grid_new = n_therm + i_cur + n_fast

    ! write energy files for debug purposes (delete)
    !open(13, file='energy_' // trim(adjustl(nuc % name)) // '_old.txt')
    !write(13, '(1e20.12)') nuc % energy
    !close(13)
    !
    !open(14, file='energy_' // trim(adjustl(nuc % name)) // '_new.txt')
    !write(14, '(1e20.12)') e_scratch(1:n_grid_new)
    !close(14)

    ! Update threshold indices. Because reactions will be accessed by a xs grid
    ! index, we offset threshold by the difference between the new xs size with
    ! n_clust values in the RRR and the old xs size with n_rrr values in the
    ! RRR.
    do i = 1, nuc % n_reaction
      rxn => nuc % reactions(i)
      if (rxn % threshold > rrr % i_high) then
        rxn % threshold = rxn % threshold + n_clust - n_rrr
      end if
    end do

    ! Update i_high index
    rrr % i_high = rrr % i_high + i_cur - n_rrr

    ! Determine offsets, which will be used to convert pointers to the energy
    ! grid to pointers to the xs grid (they are now different because one
    ! cluster may be used several times on the energy grid).
    rrr % offset_rrr = rrr %i_low - 1
    rrr % offset_fast = rrr % n_clust - (rrr % i_high - rrr % i_low + 1)

    ! Move e_scratch to energy, trimming off unused values
    n_grid_new = n_therm + i_cur + n_fast
    deallocate(nuc % energy)
    allocate(nuc % energy(n_grid_new))
    nuc % energy(1:n_grid_new) = e_scratch(1:n_grid_new)
    nuc % n_grid = n_grid_new

    ! Move c_scratch to codebook, trimming off unused values
    deallocate(rrr % codebook)
    allocate(rrr % codebook(i_cur))
    rrr % codebook(1:i_cur) = c_scratch(1:i_cur)

  end subroutine thin_energy_grid

!===============================================================================
! REORDER_CLUSTERS reorders the clusters so cluster 1 appears first in the
! energy domain, then cluster 2, etc. Takes advantage of the fact that indexing
! of codebook is the same order as indexing of the energy domain.
!===============================================================================

  subroutine reorder_clusters(codebook, n_therm, n_rrr, n_clust)

    integer, allocatable, intent(inout) :: codebook(:) ! how to do the condensation
    integer, intent(in) :: n_therm ! size of xs below the RRR
    integer, intent(in) :: n_rrr   ! size of RRR grid before clustering
    integer, intent(in) :: n_clust ! number of clusters
    integer, allocatable :: c_map(:) ! mapping from old to new cluster numbers
    logical, allocatable :: c_done(:) ! whether each cluster has been found
    integer :: offset_rrr ! offset to RRR in energy
    integer :: i ! energy index
    integer :: i_clust ! old cluster index
    integer :: i_cur ! new cluster index

    ! Determine offsets and allocate arrays
    offset_rrr = n_therm
    allocate(c_map(n_clust))
    c_map = 0
    allocate(c_done(n_clust))
    c_done = .false.

    !write (*,*) 'c_map'
    !write (*,'(30i3)') c_map

    ! Loop through to find first useage 
    i_cur = 0
    do i = 1, n_rrr
      i_clust = codebook(i)
      if (.not. c_done(i_clust)) then
        c_done(i_clust) = .true.
        i_cur = i_cur + 1
        c_map(i_clust) = i_cur
        if (i_cur == n_clust) exit
      end if
    end do

    !write (*,*) 'c_map'
    !write (*,'(30i3)') c_map

    ! Deal with empty cluster possibility
    if (i_cur /= n_clust) then
      do i_clust = 1, n_clust
        if (.not. c_done(i_clust)) then
          c_done(i_clust) = .true.
          i_cur = i_cur + 1
          c_map(i_clust) = i_cur
          if (i_cur == n_clust) exit
        end if
      end do
    end if

    !write (*,*) 'codebook'
    !write (*,'(30i3)') codebook(1:min(size(codebook),400))

    ! Map old cluster numbers to new cluster numbers
    codebook(1:n_rrr) = c_map(codebook(1:n_rrr))

    !write (*,*) 'c_map'
    !write (*,'(30i3)') c_map
    !write (*,*) 'codebook'
    !write (*,'(30i3)') codebook(1:min(size(codebook),400))

  end subroutine reorder_clusters

!===============================================================================
! GET_XS_SIZES determines and prints the current sizes of xs and energy grids.
! Unallocated arrays will print zero size.
!===============================================================================

  subroutine get_xs_sizes()
    integer(8) :: n_xs_tot   ! sum over nuclides of total cross section sizes
    integer(8) :: n_xs_main  ! sum over nuclides of tot, abs, el, (fis) sizes
    integer(8) :: n_xs_rxn   ! sum over nuclides of all reactions' sigmas
    integer(8) :: n_xs_urr   ! size of URR energy grid and probability tables
    integer(8) :: n_xs_2nd   ! size of secondary distributions for reactions
    integer(8) :: n_xs_other ! size of other xs uncounted previously
    integer(8) :: n_grid ! sum of all grid index sizes
    integer(8) :: n_e    ! sum over nuclides of the local energy grid sizes
    integer(8) :: n_ueg  ! union energy grid size
    integer(8) :: n_code ! sum of all codebook sizes
    integer(8) :: n_tot  ! total size
    integer :: s_real    ! size of real(8) in B
    integer :: s_int     ! size of int in B
    integer :: i_nuclide ! iteration index over nuclides
    integer :: i_rxn     ! iteration index over reactions
    integer :: i_2nd     ! iteration index over secondary energy distribution
    integer :: KB        ! conversion from B to KB
    integer :: MB        ! conversion from B to MB
    character(len=10) :: num ! temporary storage
    type(Nuclide), pointer :: nuc => null()
    type(Reaction), pointer :: rxn => null()
    type(RrrData), pointer :: rrr => null()

    ! Initialize to zero
    n_xs_tot = 0
    n_xs_main = 0
    n_xs_rxn = 0
    n_xs_urr = 0
    n_xs_2nd = 0
    n_xs_other = 0
    n_grid = 0
    n_e = 0
    n_code = 0
    n_ueg = 0

    ! Conversion sizes
    s_real = 8
    s_int = 4
    KB = 1024
    MB = KB * KB

    ! Loop over nuclides and accumulate sizes
    do i_nuclide = 1, n_nuclides_total
      nuc => nuclides(i_nuclide)
      ! Total
      n_xs_tot = n_xs_tot + s_real * size(nuc % total)

      ! Main
      n_xs_main = n_xs_main + s_real * &
        (size(nuc % total) + size(nuc % elastic))
      if (allocated(nuc % absorption)) then
        n_xs_main = n_xs_main + s_real * size(nuc % absorption)
      end if
      if (nuc % fissionable) then
        n_xs_main = n_xs_main + s_real * &
          (size(nuc % fission) + size(nuc % nu_fission))
      end if

      ! Reactions
      do i_rxn = 1, nuc % n_reaction
        rxn => nuc % reactions(i_rxn)
        if (allocated(rxn % sigma)) then
          n_xs_rxn = n_xs_rxn + s_real * size(rxn % sigma)
        end if
      end do

      ! Secondary
      do i_rxn = 1, nuc % n_reaction
        rxn => nuc % reactions(i_rxn)
        if (allocated(rxn % adist % energy)) then
          n_xs_2nd = n_xs_2nd + s_real * size(rxn % adist % energy)
        end if
        if (allocated(rxn % adist % type)) then
          n_xs_2nd = n_xs_2nd + s_int * size(rxn % adist % type)
        end if
        if (allocated(rxn % adist % location)) then
          n_xs_2nd = n_xs_2nd + s_int * size(rxn % adist % location)
        end if
        if (allocated(rxn % adist % data)) then
          n_xs_2nd = n_xs_2nd + s_real * size(rxn % adist % data)
        end if
        if (associated(rxn % edist)) then
          if (allocated(rxn % edist % data)) then
            n_xs_2nd = n_xs_2nd + s_real * size(rxn % edist % data)
          end if
        end if
      end do

      ! URR
      if (nuc % urr_present) then
        n_xs_urr = n_xs_urr + s_real * ( &
          size(nuc % urr_data % energy) + size(nuc % urr_data % prob, 1) * &
          size(nuc % urr_data % prob, 2) * size(nuc % urr_data % prob, 3))
      end if

      ! Other
      if (nuc % fissionable .and. allocated(nuc % nu_t_data)) then
        n_xs_other = n_xs_other + s_real * size(nuc % nu_t_data)
      end if
      if (nuc % fissionable .and. allocated(nuc % nu_p_data)) then
        n_xs_other = n_xs_other + s_real * size(nuc % nu_p_data)
      end if
      if (nuc % fissionable .and. allocated(nuc % nu_d_data)) then
        n_xs_other = n_xs_other + s_real * size(nuc % nu_d_data)
      end if
      if (nuc % fissionable .and. allocated(nuc % nu_d_precursor_data)) then
        n_xs_other = n_xs_other + s_real * size(nuc % nu_d_precursor_data)
      end if
      if (nuc % fissionable .and. associated(nuc % nu_d_edist)) then
        do i_2nd = 1, size(nuc % nu_d_edist)
          if (allocated(nuc % nu_d_edist(i_2nd) % data)) then
            n_xs_other = n_xs_other + &
              s_real * size(nuc % nu_d_edist(i_2nd) % data)
          end if
        end do
      end if

      ! Grid index
      if (allocated(nuc % grid_index)) then
        n_grid = n_grid + s_int * size(nuc % grid_index)
      end if

      ! Nuclide energy grid
      n_e = n_e + s_real * size(nuc % energy)

      ! Codebook
      if (nuc % rrr_cluster) then
        rrr => nuc % rrr_data
        if (allocated(rrr % codebook)) then
          n_code = n_code + s_int * size(rrr % codebook)
        end if
      end if
    end do

    ! Union energy grid
    if (allocated(e_grid)) then
      n_ueg = s_real * size(e_grid)
    end if

    ! Cumulative (overflow issues)
    n_tot = int(n_xs_main,8) + int(n_xs_rxn,8) + int(n_xs_urr,8) + &
        int(n_xs_2nd,8) + int(n_xs_other,8) + int(n_grid,8) + int(n_e,8) + &
        int(n_ueg,8) + int(n_code,8)

    ! Convert from B to KB
    n_xs_tot = n_xs_tot / KB
    n_xs_main = n_xs_main / KB
    n_xs_rxn = n_xs_rxn  / KB
    n_xs_urr = n_xs_urr / KB
    n_xs_2nd = n_xs_2nd / KB
    n_xs_other = n_xs_other / KB
    n_grid = n_grid / KB
    n_e = n_e / KB
    n_ueg = n_ueg / KB
    n_code = n_code / KB
    n_tot = n_tot / MB

    ! Print sizes
    write(num, '(i10)') n_xs_tot
    message = 'Size of all total xs is           ' // num // ' KB'
    call write_message(5)
    !
    write(num, '(i10)') n_xs_main
    message = 'Size of all main xs is            ' // num // ' KB'
    call write_message(5)
    !
    write(num, '(i10)') n_xs_rxn
    message = "Size of all reactions' sigmas is  " // num // ' KB'
    call write_message(5)
    !
    write(num, '(i10)') n_xs_urr
    message = 'Size of all urr data is           ' // num // ' KB'
    call write_message(5)
    !
    write(num, '(i10)') n_xs_2nd
    message = "Size of all rxns' 2ndary xs is    " // num // ' KB'
    call write_message(5)
    !
    write(num, '(i10)') n_xs_other
    message = 'Size of other fission data is     ' // num // ' KB'
    call write_message(5)
    !
    write(num, '(i10)') n_grid
    message = "Size of nuclides' index grids is  " // num // ' KB'
    call write_message(5)
    !
    write(num, '(i10)') n_e
    message = "Size of nuclides' energy grids is " // num // ' KB'
    call write_message(5)
    !
    write(num, '(i10)') n_ueg
    message = 'Size of unionized energy grid is  ' // num // ' KB'
    call write_message(5)
    !
    write(num, '(i10)') n_code
    message = 'Size of all codebooks is          ' // num // ' KB'
    call write_message(5)
    !
    write(num, '(i10)') n_tot
    message = 'Total size is                     ' // num // ' MB'
    call write_message(5)

  end subroutine get_xs_sizes

!===============================================================================
! WRITE_CLUSTERING writes the observation and cluster centers to file
!===============================================================================

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
