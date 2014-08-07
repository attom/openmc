module clustering_kmeans

  use constants
  use error,           only: fatal_error, warning
  use global
  use random_lcg,      only: prn         ! random numbers
  use string,          only: to_str

  implicit none

contains

!===============================================================================
! PERFORM_KMS performs k-means on a set of observations given an initial
! set of clusters
!===============================================================================

  subroutine perform_kms(observations, max_it, tol, clust_cen, &
    codebook, L2_err, Linf_err, err_flag)

    ! All 2-D arrays are indexed with (dimension index, feature index)
    real(8), allocatable, intent(in) :: observations(:,:) ! data to be clustered
    integer, intent(in) :: max_it                         ! max iterations
    real(8), intent(in) :: tol                            ! stopping criterion
    real(8), allocatable, intent(inout) :: clust_cen(:,:) ! init cluster locs
    integer, allocatable, intent(out) :: codebook(:)      ! cluster membership
    real(8), intent(out) :: L2_err                        ! clustering L^2 error
    real(8), intent(out) :: Linf_err                  ! clustering L^infty error
    logical, intent(out) :: err_flag                      ! true if error occurs
    real(8) :: old_L2_err, old_Linf_err                 ! old clustering errors
    real(8) :: change                                     ! change in error
    integer :: it                                         ! current iteration
    logical :: finished                                   ! clusters converged
    integer :: n_feat

    ! Initialize arrays
    n_feat = size(observations, 2)
    allocate(codebook(n_feat))

    ! Iterate the k-means algorithm until convergence is achieved
    old_L2_err = ZERO
    old_Linf_err = ZERO
    it = 0
    err_flag = .false.
    finished = .false.
    do while (.not. finished)
      call kms_compute_distances(observations, clust_cen, codebook, L2_err, &
        Linf_err)
      it = it + 1
      change = abs(L2_err - old_L2_err) / (L2_err + TINY_BIT)
      old_L2_err = L2_err
      !call print_kms(observations, clust_cen, codebook, L2_err, &
      !  it, change, tol)
      !write (*,*) it, change, L2_err, Linf_err
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
    L2_err, Linf_err)

    ! All 2-D arrays are indexed with (dimension index, feature index)
    real(8), allocatable, intent(in) :: observations(:,:) ! data to be clustered
    real(8), allocatable, intent(in) :: clust_cen(:,:)    ! cluster locs
    integer, allocatable, intent(inout) :: codebook(:)    ! cluster membership
    real(8), intent(out) :: L2_err                        ! dist data to clust
    real(8), intent(out) :: Linf_err                      ! another error metric
    integer :: n_dim                                      ! number of dimensions
    integer :: n_feat                                     ! number of features
    integer :: n_clust                                    ! number of clusters
    integer :: i_feat, i_clust
    integer :: closest_clust
    real(8) :: L2_dist, min_L2_dist, Linf_dist

    ! Sizes
    n_dim = size(observations, 1)
    n_feat = size(observations, 2)
    n_clust = size(clust_cen, 2)

    ! Compute distances from each point to each cluster. Update errors and
    ! codebook.
    L2_err = ZERO
    LInf_err = ZERO
    do i_feat = 1, n_feat
      min_L2_dist = INFINITY
      do i_clust = 1, n_clust
        L2_dist = sum((observations(:, i_feat) - clust_cen(:, i_clust))**2)
        if (L2_dist < min_L2_dist) then
          min_L2_dist = L2_dist
          closest_clust = i_clust
          Linf_dist = maxval(abs(observations(:,i_feat) - clust_cen(:, i_clust)))
        end if
      end do
      L2_err = L2_err + min_L2_dist
      Linf_err = max(Linf_err, Linf_dist)
      codebook(i_feat) = closest_clust
    end do
    L2_err = sqrt(L2_err / (n_feat * n_dim))

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
    integer :: i_feat, i_clust

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
      ! If a cluster ends up being unused, set an error
      if (clust_count(i_clust) == 0) then
        clust_count(i_clust) = 1
        clust_cen(:, i_clust) = 0
        message = "Resetting unused cluster."
        call warning()
        !err_flag = .true.
        !exit
      end if
      clust_cen(:, i_clust) = clust_cen(:, i_clust) / clust_count(i_clust)
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
      message = "Invalid number of clusters " // to_str(n_clust) // "."
      call fatal_error()
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

end module clustering_kmeans
