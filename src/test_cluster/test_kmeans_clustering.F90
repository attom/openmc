module kmeans_test_inputs

  use constants
  use global

  implicit none

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

end module kmeans_test_inputs

!===============================================================================
! MAIN runs a test of k-means
!===============================================================================

program main

  use clustering_kmeans,  only: perform_kms
  use kmeans_test_inputs

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
