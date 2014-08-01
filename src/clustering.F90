module clustering

  !use constants
  !use error,           only: fatal_error
  !use global
  !use random_lcg,      only: prn         ! random numbers

  implicit none

contains

!===============================================================================
! PERFORM_KMS performs k-means on a set of observations given an initial
! set of clusters
!===============================================================================

  subroutine perform_kms(observations, max_it, tol, clust_cen, &
    codebook, distortion, err_flag)

    ! All 2-D arrays are indexed with (dimension index, feature index)
    real(8), allocatable, intent(in) :: observations(:,:) ! data to be clustered
    integer, intent(in) :: max_it                         ! max iterations
    real(8), intent(in) :: tol                            ! stopping criterion
    real(8), allocatable, intent(inout) :: clust_cen(:,:) ! init cluster locs
    integer, allocatable, intent(out) :: codebook(:)      ! cluster membership
    real(8), intent(out) :: distortion                    ! clustering error
    logical, intent(out) :: err_flag                      ! true if error occurs
    real(8) :: old_distortion                             ! old clustering error
    real(8) :: change                                     ! change in distortion
    integer :: it                                         ! current iteration
    logical :: finished                                   ! clusters converged
    real(8) :: ZERO = 0.0_8 ! delete
    real(8) :: TINY_BIT = 1e-8_8 !delete
    integer :: n_feat

    ! Initialize arrays
    n_feat = size(observations, 2)
    allocate(codebook(n_feat))

    ! Iterate the k-means algorithm until convergence is achieved
    old_distortion = ZERO
    it = 0
    err_flag = .false.
    finished = .false.
    do while (.not. finished)
      call kms_compute_distances(observations, clust_cen, codebook, distortion)
      it = it + 1
      change = abs(distortion - old_distortion) / (distortion + TINY_BIT)
      old_distortion = distortion
      call print_kms(observations, clust_cen, codebook, distortion, &
        it, change, tol)
      call kms_update_clust_cen(observations, codebook, clust_cen, err_flag)
      if (it >= max_it .or. change <= tol .or. err_flag) then
        finished = .true.
      end if
    end do

  end subroutine perform_kms

!===============================================================================
! KMS_COMPUTE_DISTANCES compute distances from each feature to its cluster,
! compute distortion and update cluster memberships
!===============================================================================

  subroutine kms_compute_distances(observations, clust_cen, codebook, &
    distortion)

    ! All 2-D arrays are indexed with (dimension index, feature index)
    real(8), allocatable, intent(in) :: observations(:,:) ! data to be clustered
    real(8), allocatable, intent(in) :: clust_cen(:,:)    ! cluster locs
    integer, allocatable, intent(inout) :: codebook(:)    ! cluster membership
    real(8), intent(out) :: distortion                    ! dist data to clust
    integer :: n_dim                                      ! number of dimensions
    integer :: n_feat                                     ! number of features
    integer :: n_clust                                    ! number of clusters
    integer :: i_dim, i_feat, i_clust
    integer :: closest_clust
    real(8) :: dist, min_dist
    real(8) :: ZERO = 0.0_8 ! delete
    real(8) :: INFINITY = huge(0.0_8) ! delete

    ! Sizes
    n_dim = size(observations, 1)
    n_feat = size(observations, 2)
    n_clust = size(clust_cen, 2)

    ! Compute distances from each point to each cluster. Update distortion and
    ! codebook.
    distortion = ZERO
    do i_feat = 1, n_feat
      min_dist = INFINITY
      do i_clust = 1, n_clust
        dist = sum((observations(:, i_feat) - clust_cen(:, i_clust))**2)
        if (dist < min_dist) then
          min_dist = dist
          closest_clust = i_clust
        end if
      end do
      distortion = distortion + min_dist
      codebook(i_feat) = closest_clust
    end do
    distortion = sqrt(distortion)

  end subroutine kms_compute_distances


!===============================================================================
! KMS_UPDATE_CLUST_LOC make each cluster the midpoint of its members
!===============================================================================

  subroutine kms_update_clust_cen(observations, codebook, clust_cen, err_flag)

    ! All 2-D arrays are indexed with (dimension index, feature index)
    real(8), allocatable, intent(in) :: observations(:,:) ! data to be clustered
    integer, allocatable, intent(in) :: codebook(:)       ! cluster membership
    real(8), allocatable, intent(inout) :: clust_cen(:,:) ! cluster locs
    logical, intent(out) :: err_flag                      ! true if error occurs
    integer, allocatable :: clust_count(:) ! number of features in each cluster
    integer :: n_dim                                      ! number of dimensions
    integer :: n_feat                                     ! number of features
    integer :: n_clust                                    ! number of clusters
    integer :: i_dim, i_feat, i_clust
    real(8) :: ZERO = 0.0_8 ! delete
    real(8) :: INFINITY = huge(0.0_8) ! delete

    ! Sizes
    n_dim = size(observations, 1)
    n_feat = size(observations, 2)
    n_clust = size(clust_cen, 2)
    allocate(clust_count(n_clust))
    clust_count = 0

    ! Compute average location for each cluster
    err_flag = .false.
    clust_cen = 0
    do i_feat = 1, n_feat
      i_clust = codebook(i_feat)
      clust_cen(:, i_clust) = clust_cen(:, i_clust) + observations(:, i_feat)
      clust_count(i_clust) = clust_count(i_clust) + 1
    end do
    do i_clust = 1, n_clust
      clust_cen(:, i_clust) = clust_cen(:, i_clust) / clust_count(i_clust)
      ! If a cluster ends up being unused, set an error
      if (clust_count(i_clust) == 0) then
        err_flag = .true.
        exit
      end if
    end do

  end subroutine kms_update_clust_cen

!===============================================================================
! KMS_UNIFORM_INITIAL_LOC creates the initial points uniformly from 
! observations
!===============================================================================

  subroutine kms_uniform_clust_cen(observations, n_clust, clust_cen)

    ! All 2-D arrays are indexed with (dimension index, feature index)
    real(8), allocatable, intent(in) :: observations(:,:) ! data to be clustered
    integer, intent(in) :: n_clust                        ! number of clusters
    real(8), allocatable, intent(out) :: clust_cen(:,:)   ! init cluster locs
    integer :: n_dim                                      ! number of dimensions
    integer :: n_feat                                     ! number of features
    integer :: stride  ! distance between initial cluster locations

    ! Sizes
    n_dim = size(observations, 1)
    n_feat = size(observations, 2)

    ! Uniformly sample observations for the initial points
    stride = n_feat / n_clust
    if (mod(n_feat, n_clust) /= 0) then
      stride = stride + 1
    end if
    if (stride == 0) then
      !message = "Invalid number of clusters " // to_str(n_clust   ) // "."
      !call fatal_error()
      stop
    end if

    allocate(clust_cen(n_dim, n_clust))
    clust_cen = observations(:, 1:n_feat:stride)

  end subroutine kms_uniform_clust_cen

!===============================================================================
! PRINT_KMS prints debug information
!===============================================================================

  subroutine print_kms(observations, clust_cen, codebook, &
    distortion, it, change, tol)

    ! All 2-D arrays are indexed with (dimension index, feature index)
    real(8), allocatable, intent(in) :: observations(:,:) ! data to be clustered
    real(8), allocatable, intent(in) :: clust_cen(:,:)    ! cluster locs
    integer, allocatable, intent(in) :: codebook(:)       ! cluster membership
    real(8), intent(in) :: distortion                     ! dist data to clust
    integer, intent(in) :: it                             ! iteration count
    real(8), intent(in) :: change                         ! rel dif
    real(8), intent(in) :: tol                            ! tolerance
    integer :: n_dim                                      ! number of dimensions
    integer :: n_feat                                     ! number of features
    integer :: n_clust                                    ! number of clusters
    integer :: i_dim, i_feat, i_clust

    ! Sizes
    n_dim = size(observations, 1)
    n_feat = size(observations, 2)
    n_clust = size(clust_cen, 2)

    write (*,*) '------------------------------------------'
    write (*,*) "iteration", it
    write (*,*) "error", change
    write (*,*) "tolerance", tol
    write (*,*) "observations"
    do i_feat = 1, n_feat
      write(*, '(100f7.2)') (observations(i_dim, i_feat), i_dim = 1, n_dim)
    end do
    write (*,*) "cluster locations"
    do i_clust = 1, n_clust
      write(*, '(100f7.2)') (clust_cen(i_dim, i_clust), i_dim = 1, n_dim)
    end do
    write (*,*) "codebook"
    write(*, '(100i3)') (codebook(i_feat), i_feat = 1, n_feat)
    write(*,*) "distortion ", distortion

  end subroutine print_kms

!===============================================================================
! KMS_OBS_1 returns example observations for case 1
!===============================================================================

  subroutine kms_obs_1(observations, n_dim, n_feat, n_clust)

    ! All 2-D arrays are indexed with (dimension index, feature index)
    real(8), allocatable, intent(out) :: observations(:,:)! data to be clustered
    integer, intent(out) :: n_dim                         ! number of dimensions
    integer, intent(out) :: n_feat                        ! number of features
    integer, intent(out) :: n_clust                       ! number of features
    integer :: i_dim, i_feat

    ! Observation / cluster parameters
    n_dim = 1
    n_feat = 6
    n_clust = 2

    ! Initialize observations
    allocate(observations(n_dim, n_feat))
    do i_feat = 1, n_feat
      do i_dim = 1, n_dim
        observations(i_dim, i_feat) = mod(i_feat, 3)
      end do
    end do

  end subroutine kms_obs_1

!===============================================================================
! KMS_OBS_2 returns example observations for case 2
!===============================================================================

  subroutine kms_obs_2(observations, n_dim, n_feat, n_clust)

    ! All 2-D arrays are indexed with (dimension index, feature index)
    real(8), allocatable, intent(out) :: observations(:,:)! data to be clustered
    integer, intent(out) :: n_dim                         ! number of dimensions
    integer, intent(out) :: n_feat                        ! number of features
    integer, intent(out) :: n_clust                       ! number of features
    integer :: i_dim, i_feat

    ! Observation / cluster parameters
    n_dim = 3
    n_feat = 10
    n_clust = 2

    ! Initialize observations
    allocate(observations(n_dim, n_feat))
    do i_feat = 1, n_feat
      do i_dim = 1, n_dim
        observations(i_dim, i_feat) = n_feat - i_feat + (i_dim - 1.0_8) / n_dim
      end do
    end do

  end subroutine kms_obs_2

!===============================================================================
! KMS_OBS_3 returns example observations for case 3
!===============================================================================

  subroutine kms_obs_3(observations, n_dim, n_feat, n_clust)

    ! All 2-D arrays are indexed with (dimension index, feature index)
    real(8), allocatable, intent(out) :: observations(:,:)! data to be clustered
    integer, intent(out) :: n_dim                         ! number of dimensions
    integer, intent(out) :: n_feat                        ! number of features
    integer, intent(out) :: n_clust                       ! number of features
    integer :: i_dim, i_feat

    ! Observation / cluster parameters
    n_dim = 2
    n_feat = 10
    n_clust = 3

    ! Initialize observations
    allocate(observations(n_dim, n_feat))
    do i_feat = 1, n_feat
      do i_dim = 1, n_dim
        observations(i_dim, i_feat) = i_feat * i_feat * i_feat + (i_dim - 1.0_8) / n_dim
      end do
    end do

  end subroutine kms_obs_3

!===============================================================================
! KMS_OBS_4 returns example observations for case 4
!===============================================================================

  subroutine kms_obs_4(observations, n_dim, n_feat, n_clust)

    ! All 2-D arrays are indexed with (dimension index, feature index)
    real(8), allocatable, intent(out) :: observations(:,:)! data to be clustered
    integer, intent(out) :: n_dim                         ! number of dimensions
    integer, intent(out) :: n_feat                        ! number of features
    integer, intent(out) :: n_clust                       ! number of features
    integer :: i_dim, i_feat
    real(8) :: PI = 3.1415926535898_8 ! delete

    ! Observation / cluster parameters
    n_dim = 2
    n_feat = 10
    n_clust = 4

    ! Initialize observations
    allocate(observations(n_dim, n_feat))
    observations = 0
    do i_feat = 1, n_feat
      observations(1, i_feat) = cos(2 * PI * i_feat / real(n_feat))
      observations(2, i_feat) = sin(2 * PI * i_feat / real(n_feat))
    end do

  end subroutine kms_obs_4

end module clustering

!===============================================================================
! MAIN runs a test of k-means
!===============================================================================

program main

  use clustering

  implicit none
  real(8), allocatable :: observations(:,:) ! data to be clustered
  real(8), allocatable :: clust_cen(:,:)    ! init cluster locs
  integer, allocatable :: codebook(:)       ! cluster membership
  logical :: err_flag                       ! whether an error occured
  integer :: n_dim                          ! number of dimensions
  integer :: n_feat                         ! number of features
  integer :: n_clust                        ! number of clusters
  integer :: i_dim, i_feat, i_clust
  real(8) :: tol
  real(8) :: distortion
  integer :: max_it

  ! Iteration parameters
  tol = 1E-3_8
  max_it = 10

  !call kms_obs_1(observations, n_dim, n_feat, n_clust)
  call kms_obs_2(observations, n_dim, n_feat, n_clust)
  !call kms_obs_3(observations, n_dim, n_feat, n_clust)
  !call kms_obs_4(observations, n_dim, n_feat, n_clust)

  ! Get initial cluster locations and then minimize distortion by applying the
  ! k-means algorithm
  call kms_uniform_clust_cen(observations, n_clust, clust_cen)

    write (*,*) '------------------------------------------'
    write (*,*) "initially"
    write (*,*) "observations"
    do i_feat = 1, n_feat
      write(*, '(100f7.2)') (observations(i_dim, i_feat), i_dim = 1, n_dim)
    end do
    write (*,*) "cluster locations"
    do i_clust = 1, n_clust
      write(*, '(100f7.2)') (clust_cen(i_dim, i_clust), i_dim = 1, n_dim)
    end do

  call perform_kms(observations, max_it, tol, clust_cen, codebook, distortion, &
    err_flag)
  write (*,*) err_flag

end program main
