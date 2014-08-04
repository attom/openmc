module clustering

  use ace_header,        only: Nuclide
  use clustering_kmeans, only: perform_kms, kms_uniform_clust_cen
  use constants
  use error,             only: fatal_error
  use global

  implicit none

contains

!===============================================================================
! CLUSTER_ONE_NUCLIDE clusters the cross sections of one nuclide in the resolved
! resonance region
!===============================================================================

  subroutine cluster_one_nuclide(i_nuclide)

    integer :: i_nuclide ! index of nuclide to be clustered
    type(Nuclide), pointer :: nuc => null() ! pointer to nuclide
    real(8), allocatable :: scratch_xs(:) ! scratch space for one cross section

    !(observations, max_it, tol, clust_cen, codebook, distortion, err_flag)

    nuc => nuclides(i_nuclide)
    write (*,*) "Clustering nuclide ", nuc % name, " into ", n_clust_glob, &
      " clusters and ", n_group_glob, " groups."

    ! Allocate RRR data for nuclide
    nuc % rrr_cluster = .true.
    allocate(nuc % rrr_data)

  end subroutine

end module clustering
