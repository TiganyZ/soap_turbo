! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! HND X
! HND X   soap_turbo
! HND X
! HND X   soap_turbo is copyright (c) 2019-2023, Miguel A. Caro
! HND X
! HND X   soap_turbo is published and distributed under the
! HND X      Academic Software License v1.0 (ASL)
! HND X
! HND X   soap_turbo is distributed in the hope that it will be useful for non-commercial
! HND X   academic research, but WITHOUT ANY WARRANTY; without even the implied
! HND X   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! HND X   ASL for more details.
! HND X
! HND X   You should have received a copy of the ASL along with this program
! HND X   (e.g. in a LICENSE.md file); if not, you can write to the original
! HND X   licensor, Miguel Caro (mcaroba@gmail.com). The ASL is also published at
! HND X   http://github.com/gabor1/ASL
! HND X
! HND X   When using this software, please cite the following reference:
! HND X
! HND X   Miguel A. Caro. Phys. Rev. B 100, 024112 (2019)
! HND X
! HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

module soap_turbo_desc

! Kind parameters, declared locally on purpose: soap_turbo is a submodule with
! its own upstream and must not depend on the parent tree's kinds module.
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64


  use soap_turbo_radial
  use soap_turbo_angular

! private so a module that uses another does not inherit a second copy
  integer, parameter, private :: sp  = real32
  integer, parameter, private :: dp  = real64
  integer, parameter, private :: i32 = int32
  integer, parameter, private :: i64 = int64
  private :: real32, real64, int32, int64

! ---------------------------------------------------------------------------
! Persistent work buffers.
!
! get_soap is called once per batch of sites, and the arrays whose second
! dimension is n_atom_pairs are large: cnk_rad_der and friends are
! (k_max, n_max, n_atom_pairs) complex, a few hundred MB apiece for l_max = 8.
! Allocating and freeing them per call means faulting in the better part of a
! gigabyte of fresh pages on every batch, which is where the time attributed to
! "Memory allocation" actually went -- not in malloc, in the page faults.
!
! Keeping them here, grown but never shrunk, faults those pages once for the
! whole run. The peak footprint is unchanged, which is what the batching is
! there to bound. Views of the exact shape are taken with pointer rank
! remapping, so the code below is written against the same shapes as before.
! ---------------------------------------------------------------------------
  real(dp),    allocatable, save, target, private :: wr_rad(:), wr_rad_der(:)
  real(dp),    allocatable, save, target, private :: wr_soap_rad(:), wr_soap_azi(:), wr_soap_pol(:)
  complex(dp), allocatable, save, target, private :: wc_ang(:), wc_ang_rad(:), wc_ang_azi(:), wc_ang_pol(:)
  complex(dp), allocatable, save, target, private :: wc_cnk(:)

  contains

!**************************************************************************
! Grow-only reallocation of the persistent buffers above.
  subroutine ensure_real_buffer(buf, n)
    implicit none
    real(dp), allocatable, intent(inout) :: buf(:)
    integer(i64), intent(in) :: n
    if( allocated(buf) )then
      if( size(buf, kind=i64) >= n ) return
      deallocate( buf )
    end if
    allocate( buf(1:n) )
  return
  end subroutine
!**************************************************************************

!**************************************************************************
  subroutine ensure_complex_buffer(buf, n)
    implicit none
    complex(dp), allocatable, intent(inout) :: buf(:)
    integer(i64), intent(in) :: n
    if( allocated(buf) )then
      if( size(buf, kind=i64) >= n ) return
      deallocate( buf )
    end if
    allocate( buf(1:n) )
  return
  end subroutine
!**************************************************************************

!**************************************************************************
  subroutine get_soap(n_sites, n_neigh, n_species, species, species_multiplicity, n_atom_pairs, mask, rjs, &
                      thetas, phis, alpha_max, l_max, rcut_hard, rcut_soft, nf, global_scaling, atom_sigma_r, &
                      atom_sigma_r_scaling, atom_sigma_t, atom_sigma_t_scaling, &
                      amplitude_scaling, radial_enhancement, central_weight, basis, scaling_mode, do_timing, &
                      do_derivatives, compress_soap, compress_P_nonzero, compress_P_i, compress_P_j, &
                      compress_P_el, soap, soap_cart_der)

  implicit none

!-------------------
! Input variables
  real(dp), intent(in) :: rjs(:), thetas(:), phis(:)
  real(dp), intent(in) :: amplitude_scaling(:), atom_sigma_r_scaling(:), atom_sigma_t(:), atom_sigma_t_scaling(:)
  real(dp), intent(in) :: central_weight(:), atom_sigma_r(:), global_scaling(:)
  real(dp), intent(in) :: nf(:), rcut_hard(:), rcut_soft(:)
  real(dp), intent(in) :: compress_P_el(:)

  integer, intent(in) :: n_species, radial_enhancement, species(:,:), species_multiplicity(:)
  integer, intent(in) :: n_sites, n_neigh(:), l_max, n_atom_pairs, alpha_max(:), compress_P_nonzero, &
                         compress_P_i(:), compress_P_j(:)

  logical, intent(in) :: do_derivatives, do_timing, mask(:,:), compress_soap

  character(*), intent(in) :: basis, scaling_mode

! Output variables
  real(dp), intent(inout) :: soap(:,:), soap_cart_der(:,:,:)
!-------------------


!-------------------
! Internal variables
  complex(dp), pointer, contiguous :: angular_exp_coeff(:,:), cnk(:,:,:)
  complex(dp), pointer, contiguous :: angular_exp_coeff_rad_der(:,:), angular_exp_coeff_azi_der(:,:)
  complex(dp), pointer, contiguous :: angular_exp_coeff_pol_der(:,:)
  complex(dp), allocatable :: eimphi(:), prefm(:), eimphi_rad_der(:)

  real(dp), allocatable, save :: W(:,:), S(:,:), multiplicity_array(:), atom_sigma_r_prev(:), rcut_hard_prev(:)
  real(dp), pointer, contiguous :: soap_rad_der(:,:), soap_azi_der(:,:), soap_pol_der(:,:)
  real(dp), allocatable :: sqrt_dot_p(:)
  real(dp), allocatable :: W_temp(:,:), S_temp(:,:)
  real(dp), pointer, contiguous :: radial_exp_coeff(:,:), radial_exp_coeff_der(:,:)
  real(dp), allocatable :: preflm(:), plm_array(:), prefl(:), fact_array(:), prefl_rad_der(:)
  real(dp), allocatable :: this_soap(:)
  real(dp) :: amplitude, multiplicity, pi, rcut_max, fourpi_rad
! Real and imaginary parts held explicitly, so that the SOAP contractions below
! can be written as the dot products they are instead of as complex products
! whose imaginary half is then thrown away.
  real(dp) :: un_r, un_i, up_r, up_i, dn_r, dn_i, dp_r, dp_i
  real(dp) :: acc_val, acc_rad, acc_azi, acc_pol
  real(dp) :: inv_p, inv_p3, dot_r, dot_a, dot_o, sv, vr, va, vo
  real(dp) :: sin_t, cos_t, sin_f, cos_f, inv_r, cx, cy, cz, px, py, pz, ax, ay
  real(dp) :: radial_time, angular_time, coeff_time, time3, total_time, soap_time, time1, time2, compress_time, &
            memory_time, basis_time

  integer, allocatable :: i_beg(:), i_end(:)
  integer, allocatable, save :: alpha_max_prev(:)
  integer, save :: n_max_prev, l_max_prev
  integer :: k_max, n_max
  integer(i64) :: nb_pair
  integer :: i, counter, j, k, n_soap, k2, k3, n, l, m, np, counter2, n_soap_uncompressed, kbase
  logical, allocatable :: do_central(:), skip_soap_component(:,:,:), skip_soap_component_flattened(:)
  logical, allocatable, save :: skip_soap_component_flattened_prev(:)
  logical, save :: recompute_basis = .true., recompute_multiplicity_array = .true.
  character*64, save :: basis_prev = "None"
!-------------------

  if( do_timing )then
    call cpu_time(time3)
  end if

!-------------------
! Constants
  pi = dacos(-1.d0)
!-------------------



  if( do_timing )then
    call cpu_time(time1)
  end if

! Check if we need to expand the central atom
  allocate( do_central(1:n_species) )
  do i = 1, n_species
    if( central_weight(i) /= 0.d0 )then
      do_central(i) = .true.
    else
      do_central(i) = .false.
    end if
  end do


  rcut_max = 0.d0
  do i = 1, n_species
    if( rcut_hard(i) > rcut_max )then
      rcut_max = rcut_hard(i)
    end if
  end do


! This assigns begin and end indices to the components of the radial basis
  allocate( i_beg(1:n_species) )
  allocate( i_end(1:n_species) )
  i_beg(1) = 1
  i_end(1) = alpha_max(1)
  do i = 2, n_species
    i_beg(i) = i_end(i-1) + 1
    i_end(i) = i_beg(i) + alpha_max(i) - 1
  end do


! Do some array allocation for the SOAP expansion part
  k_max = 1 + l_max*(l_max+1)/2 + l_max
  allocate( preflm(1:k_max) )
  allocate( plm_array(1:k_max) )
  allocate( eimphi(1:k_max) )
  allocate( prefl(0:l_max) )
  allocate( prefm(0:l_max) )
  allocate( fact_array(1:l_max) )
  call get_preflm(preflm, l_max)
  if( do_derivatives )then
    allocate( prefl_rad_der(0:l_max) )
    allocate( eimphi_rad_der(1:k_max) )
  end if


  if( do_timing )then
    call cpu_time(time2)
    memory_time = time2 - time1
    time1 = time2
  end if


! This is to build the radial basis
  n_max = 0
  do i = 1, n_species
    n_max = n_max + alpha_max(i)
  end do
! In many cases the basis does not need to be recomputed, however if the basis definition
! changes, the matrices that relate the primitive (polynomial) basis to the orthonormal
! basis need to be recomputed. For the purely polynomial bases (poly3, poly3tab and
! poly3operator), only the number of basis functions matters. For the Gaussian augmented
! basis (poly3gauss), the basis depends also on the atom_sigma_r and rcut_hard values
!
! Trivial
  if( trim(adjustl(basis)) /= trim(adjustl(basis_prev)) )then
    recompute_basis = .true.
    basis_prev = trim(adjustl(basis))
  end if
! Check for changes in number of basis functions: recompute basis for all basis types
  if( .not. allocated(alpha_max_prev) )then
    allocate( alpha_max_prev(1:size(alpha_max)) )
    recompute_basis = .true.
  else if( size(alpha_max) /= size(alpha_max_prev) )then
    deallocate( alpha_max_prev )
    allocate( alpha_max_prev(1:size(alpha_max)) )
    recompute_basis = .true.
  else if( any( alpha_max /= alpha_max_prev ) )then
    recompute_basis = .true.
  end if
! Check for changes in atom_sigma_r: recompute basis only for poly3gauss
  if( .not. allocated(atom_sigma_r_prev) )then
    allocate( atom_sigma_r_prev(1:size(atom_sigma_r)) )
    if( basis == "poly3gauss" )then
      recompute_basis = .true.
    end if
  else if( size(atom_sigma_r) /= size(atom_sigma_r_prev) )then
    deallocate( atom_sigma_r_prev )
    allocate( atom_sigma_r_prev(1:size(atom_sigma_r)) )
    if( basis == "poly3gauss" )then
      recompute_basis = .true.
    end if
  else if( any( atom_sigma_r /= atom_sigma_r_prev ) )then
    if( basis == "poly3gauss" )then
      recompute_basis = .true.
    end if
  end if
! Check for changes in rcut_hard: recompute basis only for poly3gauss
  if( .not. allocated(rcut_hard_prev) )then
    allocate( rcut_hard_prev(1:size(rcut_hard)) )
    if( basis == "poly3gauss" )then
      recompute_basis = .true.
    end if
  else if( size(rcut_hard) /= size(rcut_hard_prev) )then
    deallocate( rcut_hard_prev )
    allocate( rcut_hard_prev(1:size(rcut_hard)) )
    if( basis == "poly3gauss" )then
      recompute_basis = .true.
    end if
  else if( any( rcut_hard /= rcut_hard_prev ) )then
    if( basis == "poly3gauss" )then
      recompute_basis = .true.
    end if
  end if
! Store for next check
  alpha_max_prev = alpha_max
  atom_sigma_r_prev = atom_sigma_r
  rcut_hard_prev = rcut_hard
  n_max_prev = n_max
! Here we compute the basis transformation matrices if they need to be recomputed
  if( recompute_basis )then
    if( allocated(W) .or. allocated(S) )then
      deallocate(W, S)
    end if
    allocate( W(1:n_max, 1:n_max) )
    allocate( S(1:n_max, 1:n_max) )
    W = 0.d0
    S = 0.d0
!   This is done per species. Each loop iteration modifies the slice of the W and S matrices that
!   corresponds to the species in question. The radial basis functions for species A are always
!   assumed orthogonal to the basis functions for species B. W and S are therefore block diagonal.
    do i = 1, n_species
!     We pass these temp arrays with the right size because the Lapack/Blas routines internally fail
!     if the memory is not contiguous
      allocate( S_temp(1:alpha_max(i), 1:alpha_max(i)) )
      allocate( W_temp(1:alpha_max(i), 1:alpha_max(i)) )
      S_temp = 0.d0
      W_temp = 0.d0
      if( basis == "poly3gauss" )then
        call get_orthonormalization_matrix_poly3gauss(alpha_max(i), atom_sigma_r(i), rcut_hard(i), S_temp, W_temp)
      else if( basis == "poly3" )then
        call get_orthonormalization_matrix_poly3(alpha_max(i), S_temp, W_temp)
      else if( basis == "poly3tab" )then
        call get_orthonormalization_matrix_poly3_tabulated(alpha_max(i), S_temp, W_temp)
      else if( basis == "poly3operator" )then
        call get_orthonormalization_matrix_poly3_tabulated(alpha_max(i), S_temp, W_temp)
      end if
      S(i_beg(i):i_end(i), i_beg(i):i_end(i)) = S_temp
      W(i_beg(i):i_end(i), i_beg(i):i_end(i)) = W_temp
      deallocate( S_temp, W_temp )
    end do
  end if
  recompute_basis = .false.

  if( do_timing )then
    call cpu_time(time2)
    basis_time = time2 - time1
  end if

! Uncompressed SOAP dimensions
  n_soap_uncompressed = n_max*(n_max+1)/2 * (l_max+1)

! This is for the expansion coefficients and the soap vectors
! radial_exp_coeff and angular_exp_coeff are not cleared here: the radial
! routines clear their own species block and the angular routine writes every
! column it does not clear itself, so a sweep here would only be a second pass
! over the same memory. cnk is accumulated over neighbours, so it does need it.
  nb_pair = int(n_atom_pairs, i64)
  call ensure_real_buffer   ( wr_rad, int(n_max,i64)*nb_pair )
  call ensure_complex_buffer( wc_ang, int(k_max,i64)*nb_pair )
  call ensure_complex_buffer( wc_cnk, int(k_max,i64)*int(n_max,i64)*int(n_sites,i64) )
  radial_exp_coeff (1:n_max, 1:n_atom_pairs)          => wr_rad(1:int(n_max,i64)*nb_pair)
  angular_exp_coeff(1:k_max, 1:n_atom_pairs)          => wc_ang(1:int(k_max,i64)*nb_pair)
  cnk              (1:k_max, 1:n_max, 1:n_sites)      => wc_cnk(1:int(k_max,i64)*int(n_max,i64)*int(n_sites,i64))
  cnk = 0.d0
! Handle SOAP compression here
  allocate( skip_soap_component_flattened(1:n_soap_uncompressed) )
  skip_soap_component_flattened = .false.
  if( compress_soap )then
    if( do_timing )then
      call cpu_time(time1)
    end if
!   New way to compress with a linear transformation
    n_soap = maxval(compress_P_i)
    skip_soap_component_flattened = .true.
    do i = 1, compress_P_nonzero
      skip_soap_component_flattened(compress_P_j(i)) = .false.
    end do

    if( do_timing )then
      call cpu_time(time2)
      compress_time = time2 - time1
    end if
  else
    n_soap = n_soap_uncompressed
  end if

  if( do_timing )then
    call cpu_time(time1)
  end if

  allocate( sqrt_dot_p(1:n_sites) )
  sqrt_dot_p = 0.d0
  if( do_derivatives )then
    call ensure_real_buffer   ( wr_rad_der,  int(n_max,i64)*nb_pair )
    call ensure_complex_buffer( wc_ang_rad,  int(k_max,i64)*nb_pair )
    call ensure_complex_buffer( wc_ang_azi,  int(k_max,i64)*nb_pair )
    call ensure_complex_buffer( wc_ang_pol,  int(k_max,i64)*nb_pair )
    call ensure_real_buffer   ( wr_soap_rad, int(n_soap,i64)*nb_pair )
    call ensure_real_buffer   ( wr_soap_azi, int(n_soap,i64)*nb_pair )
    call ensure_real_buffer   ( wr_soap_pol, int(n_soap,i64)*nb_pair )
    radial_exp_coeff_der     (1:n_max, 1:n_atom_pairs)        => wr_rad_der (1:int(n_max,i64)*nb_pair)
    angular_exp_coeff_rad_der(1:k_max, 1:n_atom_pairs)        => wc_ang_rad (1:int(k_max,i64)*nb_pair)
    angular_exp_coeff_azi_der(1:k_max, 1:n_atom_pairs)        => wc_ang_azi (1:int(k_max,i64)*nb_pair)
    angular_exp_coeff_pol_der(1:k_max, 1:n_atom_pairs)        => wc_ang_pol (1:int(k_max,i64)*nb_pair)
    soap_rad_der(1:n_soap, 1:n_atom_pairs)                    => wr_soap_rad(1:int(n_soap,i64)*nb_pair)
    soap_azi_der(1:n_soap, 1:n_atom_pairs)                    => wr_soap_azi(1:int(n_soap,i64)*nb_pair)
    soap_pol_der(1:n_soap, 1:n_atom_pairs)                    => wr_soap_pol(1:int(n_soap,i64)*nb_pair)
!   Nothing is cleared here either. The radial and angular routines clear what
!   they do not write; get_derivatives clears the cnk_*_der slab of any pair it
!   skips; and every column of soap_*_der is either assigned outright or, in
!   the compressed case, cleared just before it is accumulated into. For
!   l_max = 8 and a 20k-pair batch these ten sweeps came to well over a
!   gigabyte of stores per call, all of it overwritten immediately after.
  else
!   We need this dummy variable defined here. Note the decreased range for the second index
    call ensure_real_buffer( wr_rad_der, int(n_max,i64) )
    radial_exp_coeff_der(1:n_max, 1:1) => wr_rad_der(1:int(n_max,i64))
  end if


  if( do_timing )then
    call cpu_time(time2)
    memory_time = memory_time + time2 - time1
    time1 = time2
  end if


  do i = 1, n_species
    if( basis == "poly3gauss" )then
      call get_radial_expansion_coefficients_poly3gauss(n_sites, n_neigh, rjs, alpha_max(i), rcut_soft(i), &
                                                        rcut_hard(i), atom_sigma_r(i), atom_sigma_r_scaling(i), &
                                                        amplitude_scaling(i), nf(i), W(i_beg(i):i_end(i),i_beg(i):i_end(i)), &
                                                        scaling_mode, mask(:,i), radial_enhancement, do_derivatives, &
                                                        radial_exp_coeff(i_beg(i):i_end(i), :), &
                                                        radial_exp_coeff_der(i_beg(i):i_end(i), :) )
      radial_exp_coeff(i_beg(i):i_end(i), :) = radial_exp_coeff(i_beg(i):i_end(i), :) * global_scaling(i)
      radial_exp_coeff_der(i_beg(i):i_end(i), :) = radial_exp_coeff_der(i_beg(i):i_end(i), :) * global_scaling(i)
    else if( basis == "poly3" )then
      call get_radial_expansion_coefficients_poly3(n_sites, n_neigh, rjs, alpha_max(i), rcut_soft(i), &
                                                   rcut_hard(i), atom_sigma_r(i), atom_sigma_r_scaling(i), &
                                                   amplitude_scaling(i), nf(i), W(i_beg(i):i_end(i),i_beg(i):i_end(i)), &
                                                   scaling_mode, mask(:,i), radial_enhancement, do_derivatives, &
                                                   do_central(i), central_weight(i), &
                                                   radial_exp_coeff(i_beg(i):i_end(i), :), &
                                                   radial_exp_coeff_der(i_beg(i):i_end(i), :) )
      radial_exp_coeff(i_beg(i):i_end(i), :) = radial_exp_coeff(i_beg(i):i_end(i), :) * global_scaling(i)
      radial_exp_coeff_der(i_beg(i):i_end(i), :) = radial_exp_coeff_der(i_beg(i):i_end(i), :) * global_scaling(i)
    else if( basis == "poly3operator" )then
      call get_radial_expansion_coefficients_poly3operator(n_sites, n_neigh, rjs, alpha_max(i), rcut_soft(i), &
                                                           rcut_hard(i), atom_sigma_r(i), atom_sigma_r_scaling(i), &
                                                           amplitude_scaling(i), W(i_beg(i):i_end(i),i_beg(i):i_end(i)), &
                                                           scaling_mode, mask(:,i), radial_enhancement, do_derivatives, &
                                                           .true.,  central_weight(i), radial_exp_coeff(i_beg(i):i_end(i), :), &
                                                           radial_exp_coeff_der(i_beg(i):i_end(i), :))
      radial_exp_coeff(i_beg(i):i_end(i), :) = radial_exp_coeff(i_beg(i):i_end(i), :) * global_scaling(i)
      radial_exp_coeff_der(i_beg(i):i_end(i), :) = radial_exp_coeff_der(i_beg(i):i_end(i), :) * global_scaling(i)
    end if
  end do



  if( do_timing )then
    call cpu_time(time2)
    radial_time = time2 - time1
    time1 = time2
  end if




! For the angular expansion the masking works differently, since we do not have a species-augmented basis as in the
! radial expansion part.
  call get_angular_expansion_coefficients(n_sites, n_neigh, thetas, phis, rjs, atom_sigma_t, atom_sigma_t_scaling, &
                                          rcut_max, l_max, eimphi, preflm, plm_array, prefl, prefm, &
                                          fact_array, mask, n_species, eimphi_rad_der, &
                                          do_derivatives, prefl_rad_der, angular_exp_coeff, angular_exp_coeff_rad_der, &
                                          angular_exp_coeff_azi_der, angular_exp_coeff_pol_der )
  if( do_timing )then
    call cpu_time(time2)
    angular_time = time2 - time1
  end if
!  write(*,"(f8.3, 1X, A)") time2-time1, "seconds"






! For debugging (only gfortran)
  if( .false. )then
    write(*,*) "# Gaussian centered at ", rjs(4)
    write(*,"(A)",advance="no") "rho(x) = "
    do i = 1, alpha_max(1)
      write(*,"(A,I0,A,E16.8,A)",advance="no") "p", i, "(x) *", radial_exp_coeff(i,4), "+"
    end do
    write(*,*) "0."
  end if





!  write(*,*) "Obtaining full expansion coefficients..."
!  write(*,*)
  if( do_timing )then
    call cpu_time(time1)
  end if
  k2 = 0
  do i = 1, n_sites
    do j = 1, n_neigh(i)
      k2 = k2 + 1
!     The l and m loops that used to be here just enumerate k = 1, k_max in
!     order, since k = 1 + l*(l+1)/2 + m. Walking k directly makes the inner
!     loop contiguous in both cnk and angular_exp_coeff.
      do n = 1, n_max
!       It is messy with the prefactor in spherical harmonics but we need to be sure because of the central atom below
        fourpi_rad = 4.d0*pi * radial_exp_coeff(n, k2)
        do k = 1, k_max
          cnk(k, n, i) = cnk(k, n, i) + fourpi_rad * angular_exp_coeff(k, k2)
        end do
      end do
    end do
    do k = 1, species_multiplicity(i)
      j = species(k, i)
      if( basis == "poly3gauss" .and. central_weight(j) /= 0.d0 )then
        if( radial_enhancement == 1 )then
          amplitude = dsqrt(2.d0/pi) * atom_sigma_r(j) / rcut_hard(j)
        else if( radial_enhancement == 2 )then
          amplitude = atom_sigma_r(j)**2 / rcut_hard(j)**2
        else
          amplitude = 1.d0
        end if
        cnk(1, i_beg(j):i_end(j), i) = cnk(1, i_beg(j):i_end(j), i) + &
                                       amplitude * central_weight(j) * dsqrt(4.d0*pi) * pi**0.25d0 * &
                                       dsqrt(atom_sigma_r(j) / 2.d0) * &
                                       rcut_hard(j)**3 / atom_sigma_t(j)**2 / atom_sigma_r(j) * &
                                       matmul(W(i_beg(j):i_end(j),i_beg(j):i_end(j)), S(i_beg(j):i_end(j), i_end(j)) )
      end if
    end do
  end do

! cnk_rad_der / cnk_azi_der / cnk_pol_der used to be built here. They are not
! built at all now: they factorise into the radial and angular expansions, and
! soap_derivatives_from_expansion contracts those directly. See the comment
! there.
  if( do_timing )then
    call cpu_time(time2)
    coeff_time = time2 - time1
  end if







!  write(*,*) "Building SOAP vectors..."

  if( do_timing )then
    call cpu_time(time1)
  end if

! Create the multiplicity array
  if( any( alpha_max_prev /= alpha_max ) .or. l_max_prev /= l_max )then
    recompute_multiplicity_array = .true.
  end if
  l_max_prev = l_max
  if( allocated( skip_soap_component_flattened_prev ) )then
    if( size(skip_soap_component_flattened_prev) /= size(skip_soap_component_flattened) )then
      recompute_multiplicity_array = .true.
      deallocate( skip_soap_component_flattened_prev )
      allocate( skip_soap_component_flattened_prev(1:size(skip_soap_component_flattened)) )
      skip_soap_component_flattened_prev = skip_soap_component_flattened
    else if( any(skip_soap_component_flattened_prev .neqv. skip_soap_component_flattened) )then
      recompute_multiplicity_array = .true.
      skip_soap_component_flattened_prev = skip_soap_component_flattened
    end if
  else
    allocate( skip_soap_component_flattened_prev(1:size(skip_soap_component_flattened)) )
    skip_soap_component_flattened_prev = skip_soap_component_flattened
  end if
     
  if( recompute_multiplicity_array )then
    counter = 0
    counter2 = 0
    do n = 1, n_max
      do np = n, n_max
        do l = 0, l_max
          counter = counter + 1
          if( skip_soap_component_flattened(counter) )cycle
          do m = 0, l
            counter2 = counter2 + 1
          end do
        end do
      end do
    end do
    if( allocated(multiplicity_array) )then
      deallocate(multiplicity_array)
    end if
    allocate( multiplicity_array(1:counter2) )
    counter = 0
    counter2 = 0
    do n = 1, n_max
      do np = n, n_max
        do l = 0, l_max
          counter = counter + 1
          if( skip_soap_component_flattened(counter) )cycle
          do m = 0, l
            counter2 = counter2 + 1
            multiplicity = 1.d0
            if( n /= np )then
              multiplicity = multiplicity * dsqrt(2.d0)
            end if
            if( m > 0 )then
              multiplicity = multiplicity * 2.d0
            end if
            multiplicity_array(counter2) = multiplicity
          end do
        end do
      end do
    end do
  end if
  recompute_multiplicity_array = .false.


  allocate( this_soap(1:n_soap_uncompressed) )
! The skipped components are the only ones the loop below does not assign, and
! they are never read either -- the compressed path reads this_soap only at the
! compress_P_j indices, which are exactly the components that are not skipped.
! So they are cleared once here instead of the whole vector being cleared per
! site.
  do counter = 1, n_soap_uncompressed
    if( skip_soap_component_flattened(counter) ) this_soap(counter) = 0.d0
  end do
  do i = 1, n_sites
    counter = 0
    counter2 = 0
    do n = 1, n_max
      do np = n, n_max
        do l = 0, l_max
          counter = counter+1
          if( skip_soap_component_flattened(counter) )cycle
          kbase = 1 + l*(l+1)/2
          acc_val = 0.d0
          do m = 0, l
            k = kbase + m
            counter2 = counter2 + 1
!            multiplicity = 1.d0
!            if( n /= np )then
!              multiplicity = multiplicity * dsqrt(2.d0)
!            end if
!            if( m > 0 )then
!              multiplicity = multiplicity * 2.d0
!            end if
!           real(z1 * conjg(z2)) is Re(z1)Re(z2) + Im(z1)Im(z2); forming the
!           whole complex product and dropping its imaginary part costs twice
!           the multiplies.
            un_r = dble(cnk(k, n, i));  un_i = aimag(cnk(k, n, i))
            up_r = dble(cnk(k, np, i)); up_i = aimag(cnk(k, np, i))
            acc_val = acc_val + multiplicity_array(counter2) * ( un_r*up_r + un_i*up_i )
          end do
          this_soap(counter) = acc_val
        end do
      end do
    end do
!   Linear transformation (it assumes sparse P; it may benefit from using a sparse linear algebra routine
!   for large dimensions and BLAS routines if P is dense)
    if( compress_soap )then
      do k = 1, compress_P_nonzero
        soap(compress_P_i(k), i) = soap(compress_P_i(k), i) + compress_P_el(k) * this_soap(compress_P_j(k))
      end do
    else
      soap(:,i) = this_soap(:)
    end if
    sqrt_dot_p(i) = dsqrt(dot_product(soap(1:n_soap, i), soap(1:n_soap, i)))
!   This is to avoid NaNs when the SOAP sphere is empty
    if( sqrt_dot_p(i) < 1.d-5 )then
      sqrt_dot_p(i) = 1.d0
    end if
  end do
  deallocate( this_soap )


  if( do_derivatives )then
!   Derivatives of the SOAP descriptor in spherical coordinates
!****************************
! Uncomment for detailed timing check

! call cpu_time(time1)
!****************************
    call soap_derivatives_from_expansion(n_sites, n_neigh, n_max, l_max, n_soap, &
                                        n_soap_uncompressed, cnk, angular_exp_coeff, &
                                        angular_exp_coeff_rad_der, angular_exp_coeff_azi_der, &
                                        angular_exp_coeff_pol_der, radial_exp_coeff, &
                                        radial_exp_coeff_der, skip_soap_component_flattened, &
                                        compress_soap, compress_P_nonzero, compress_P_i, &
                                        compress_P_j, compress_P_el, soap_rad_der, soap_azi_der, &
                                        soap_pol_der)
!   Normalise and rotate to Cartesian. The normalised spherical components are
!   consumed by the rotation in the same iteration, so they are kept in
!   registers rather than written back to soap_*_der and read again; the two
!   angles are evaluated once per pair instead of once per component
!   expression; and the 1/sqrt_dot_p factors are reciprocals rather than
!   n_soap divisions apiece.
    k2 = 0
    k3 = 0
    do i = 1, n_sites
      inv_p  = 1.d0 / sqrt_dot_p(i)
      inv_p3 = 1.d0 / sqrt_dot_p(i)**3
      do j = 1, n_neigh(i)
        k2 = k2 + 1
        if( j == 1 )then
          k3 = k2
          cycle
        end if

        dot_r = 0.d0
        dot_a = 0.d0
        dot_o = 0.d0
        do k = 1, n_soap
          sv = soap(k, i)
          dot_r = dot_r + sv*soap_rad_der(k, k2)
          dot_a = dot_a + sv*soap_azi_der(k, k2)
          dot_o = dot_o + sv*soap_pol_der(k, k2)
        end do

        sin_t = dsin(thetas(k2))
        cos_t = dcos(thetas(k2))
        sin_f = dsin(phis(k2))
        cos_f = dcos(phis(k2))
        inv_r = 1.d0 / rjs(k2)
        cx =  sin_t*cos_f
        cy =  sin_t*sin_f
        cz =  cos_t
        px = -cos_t*cos_f*inv_r
        py = -cos_t*sin_f*inv_r
        pz =  sin_t*inv_r
        ax = -sin_f*inv_r
        ay =  cos_f*inv_r

!       MAKE SURE THAT THIS IS CORRECT FOR THE CENTRAL ATOM DERIVATIVES
        do k = 1, n_soap
          sv = soap(k, i)
          vr = soap_rad_der(k, k2)*inv_p - sv*inv_p3*dot_r
          va = soap_azi_der(k, k2)*inv_p - sv*inv_p3*dot_a
          vo = soap_pol_der(k, k2)*inv_p - sv*inv_p3*dot_o
          soap_cart_der(1, k, k2) = cx*vr + px*vo + ax*va
          soap_cart_der(2, k, k2) = cy*vr + py*vo + ay*va
          soap_cart_der(3, k, k2) = cz*vr + pz*vo
        end do
        do k = 1, n_soap
          soap_cart_der(1, k, k3) = soap_cart_der(1, k, k3) - soap_cart_der(1, k, k2)
          soap_cart_der(2, k, k3) = soap_cart_der(2, k, k3) - soap_cart_der(2, k, k2)
          soap_cart_der(3, k, k3) = soap_cart_der(3, k, k3) - soap_cart_der(3, k, k2)
        end do
      end do
    end do
!****************************
! Uncomment for detailed timing check
!
! call cpu_time(time2)
! write(*,*) time2-time1
!****************************
  end if

! Now we normalize the soap vectors:
  do i = 1, n_sites
    soap(1:n_soap, i) = soap(1:n_soap, i) / sqrt_dot_p(i)
  end do


! This is for debugging
 if( .false. )then
!   open(unit=10, file="soap_desc.dat", status="unknown", access="append")
   open(unit=10, file="soap_desc.dat", status="unknown")
   write(10, *) soap(1:n_soap, 1)
   close(10)
   open(unit=10, file="soap_rad_der.dat", status="unknown", access="append")
   write(10, *) soap_rad_der(1:n_soap, 2)
   close(10)
   open(unit=10, file="soap_azi_der.dat", status="unknown", access="append")
   write(10, *) soap_azi_der(1:n_soap, 2)
   close(10)
   open(unit=10, file="soap_pol_der.dat", status="unknown", access="append")
   write(10, *) soap_pol_der(1:n_soap, 2)
   close(10)
   open(unit=10, file="soap_cart_der.dat", status="unknown")
   do i = 1, 5
     write(10, *) soap_cart_der(1, 1:n_soap, i)
     write(10, *) soap_cart_der(2, 1:n_soap, i)
     write(10, *) soap_cart_der(3, 1:n_soap, i)
     write(10, *)
   end do
   close(10)
!   open(unit=10, file="norm.dat", status="unknown", access="append")
!   write(10, *) sqrt_dot_p(1)
!   close(10)

  end if

!  call cpu_time(time2)
!  write(*,"(f8.3, 1X, A)") time2-time1, "seconds"


! The arrays sized by n_atom_pairs live in the persistent buffers above and are
! deliberately not freed here; the pointers into them are re-associated on the
! next call.
  deallocate( eimphi, preflm, plm_array, prefl, prefm, fact_array, &
              i_beg, i_end, do_central, sqrt_dot_p )
  deallocate( skip_soap_component_flattened )
  if( do_derivatives )then
    deallocate( eimphi_rad_der, prefl_rad_der )
  end if

  if( do_timing )then
    call cpu_time(time2)
    soap_time = time2-time1
    total_time = time2 - time3
    write(*,*)'                                       |'
    write(*,*)'SOAP timings:                          |'
    write(*,*)'                                       |'
    write(*,'(A, F8.3, A)') '  *) Radial expansion: ', radial_time, ' seconds |'
    write(*,'(A, F6.3, A)') '  *) Radial basis build: ', basis_time, ' seconds |'
    write(*,'(A, F7.3, A)') '  *) Angular expansion: ', angular_time, ' seconds |'
    write(*,'(A, F7.3, A)') '  *) Expansion coeffs.: ', coeff_time, ' seconds |'
    write(*,'(A, F7.3, A)') '  *) SOAP vector build: ', soap_time, ' seconds |'
    if( compress_soap )then
      write(*,'(A, F8.3, A)') '  *) SOAP compression: ', compress_time, ' seconds |'
    end if
    write(*,'(A, F7.3, A)') '  *) Memory allocation: ', memory_time, ' seconds |'
    write(*,'(A, F19.3, A)') '  *) Total: ', total_time, ' seconds |'
    write(*,*)'                                       |'
    write(*,*)'.......................................|'
  end if

  end subroutine get_soap


!**************************************************************************
!
! Spherical derivatives of the (uncompressed, unnormalised) SOAP vector,
! straight from the radial and angular expansions.
!
! The old route went through cnk_rad_der / cnk_azi_der / cnk_pol_der, three
! (k_max, n_max, n_atom_pairs) complex arrays. For l_max = 8 and n_max = 16
! that is 46 kB per neighbour pair to write and then read back, a few hundred
! MB per batch, and it is the reason "Memory allocation" and "SOAP vector
! build" between them were 80% of get_soap.
!
! They are not needed, because cnk_alpha_der factorises. With
!   A(k)   = angular_exp_coeff(k,k2),      A_r(k) = angular_exp_coeff_rad_der(k,k2)
!   r(n)   = radial_exp_coeff(n,k2),       r'(n) = radial_exp_coeff_der(n,k2)
! the old array was
!   cnk_rad_der(k,n,k2) = 4pi [ r'(n) A(k) + r(n) A_r(k) ]
! which is separable in n and k, and the SOAP derivative contracts it against
! cnk over k. Writing the multiplicity as base(n,np)*w(m), with
! base = sqrt(2) for n /= np and w = 2 for m > 0, and using
! Re(z conjg(w)) = Re(w conjg(z)), the contraction collapses to
!
!   P (n,l) = sum_m w(m) Re( A(k)     conjg(cnk(k,n,i)) )
!   Q (n,l) = sum_m w(m) Re( A_r(k)   conjg(cnk(k,n,i)) )
!   Qa(n,l) = sum_m w(m) Re( A_azi(k) conjg(cnk(k,n,i)) )
!   Qp(n,l) = sum_m w(m) Re( A_pol(k) conjg(cnk(k,n,i)) )
!
!   d_rad(n,np,l) = 4pi base [ r'(n) P(np,l) + r(n) Q(np,l)
!                            + r'(np) P(n,l) + r(np) Q(n,l) ]
!   d_azi(n,np,l) = 4pi base [ r(n) Qa(np,l) + r(np) Qa(n,l) ]
!   d_pol(n,np,l) = 4pi base [ r(n) Qp(np,l) + r(np) Qp(n,l) ]
!
! so the per-pair traffic drops from 46 kB to about 3 kB, and the n_max*(n_max+1)/2
! contractions over k become n_max of them.
!
! Pairs beyond the cutoff need no special case: the angular routine leaves
! their A columns zero, so every P and Q comes out zero, exactly as the
! explicit cnk_*_der = 0 did.
!
  subroutine soap_derivatives_from_expansion(n_sites, n_neigh, n_max, l_max, n_soap, &
                                             n_soap_uncompressed, cnk, angular_exp_coeff, &
                                             angular_exp_coeff_rad_der, angular_exp_coeff_azi_der, &
                                             angular_exp_coeff_pol_der, radial_exp_coeff, &
                                             radial_exp_coeff_der, skip, compress_soap, &
                                             compress_P_nonzero, compress_P_i, compress_P_j, &
                                             compress_P_el, soap_rad_der, soap_azi_der, soap_pol_der)

    implicit none

    integer, intent(in) :: n_sites, n_neigh(:), n_max, l_max, n_soap, n_soap_uncompressed
    integer, intent(in) :: compress_P_nonzero, compress_P_i(:), compress_P_j(:)
    complex(dp), intent(in), contiguous :: cnk(:,:,:), angular_exp_coeff(:,:)
    complex(dp), intent(in), contiguous :: angular_exp_coeff_rad_der(:,:), angular_exp_coeff_azi_der(:,:)
    complex(dp), intent(in), contiguous :: angular_exp_coeff_pol_der(:,:)
    real(dp), intent(in), contiguous :: radial_exp_coeff(:,:), radial_exp_coeff_der(:,:)
    real(dp), intent(in) :: compress_P_el(:)
    logical, intent(in) :: skip(:), compress_soap
    real(dp), intent(inout), contiguous :: soap_rad_der(:,:), soap_azi_der(:,:), soap_pol_der(:,:)

    integer :: i, j, k, k2, kb, l, m, n, np, counter, k_max
    real(dp) :: pi, fourpi, base, sq2, w
    real(dp) :: cr, ci, acc_p, acc_q, acc_a, acc_o
    real(dp) :: rn, rdn, rnp, rdnp, pre

    real(dp) :: P(n_max, 0:l_max), Q(n_max, 0:l_max), Qa(n_max, 0:l_max), Qo(n_max, 0:l_max)
    real(dp), allocatable :: this_rad(:), this_azi(:), this_pol(:)

    pi = dacos(-1.d0)
    fourpi = 4.d0*pi
    sq2 = dsqrt(2.d0)
    k_max = 1 + l_max*(l_max+1)/2 + l_max

    allocate( this_rad(1:n_soap_uncompressed) )
    allocate( this_azi(1:n_soap_uncompressed) )
    allocate( this_pol(1:n_soap_uncompressed) )
!   Skipped components are never assigned below and never read; clear them once.
    do counter = 1, n_soap_uncompressed
      if( skip(counter) )then
        this_rad(counter) = 0.d0
        this_azi(counter) = 0.d0
        this_pol(counter) = 0.d0
      end if
    end do

    k2 = 0
    do i = 1, n_sites
      do j = 1, n_neigh(i)
        k2 = k2 + 1
!       j = 1 is the central atom. Its Cartesian column is built purely by
!       accumulation from its neighbours, so its spherical derivative is never
!       read and is not computed here.
        if( j == 1 ) cycle

!       ---- the four contractions over k, one per (n,l)
        do n = 1, n_max
          do l = 0, l_max
            kb = 1 + l*(l+1)/2
            acc_p = 0.d0
            acc_q = 0.d0
            acc_a = 0.d0
            acc_o = 0.d0
            do m = 0, l
              k = kb + m
              w = 2.d0
              if( m == 0 ) w = 1.d0
              cr = dble(cnk(k,n,i))
              ci = aimag(cnk(k,n,i))
              acc_p = acc_p + w*( dble(angular_exp_coeff(k,k2))*cr + aimag(angular_exp_coeff(k,k2))*ci )
              acc_q = acc_q + w*( dble(angular_exp_coeff_rad_der(k,k2))*cr &
                                + aimag(angular_exp_coeff_rad_der(k,k2))*ci )
              acc_a = acc_a + w*( dble(angular_exp_coeff_azi_der(k,k2))*cr &
                                + aimag(angular_exp_coeff_azi_der(k,k2))*ci )
              acc_o = acc_o + w*( dble(angular_exp_coeff_pol_der(k,k2))*cr &
                                + aimag(angular_exp_coeff_pol_der(k,k2))*ci )
            end do
            P (n,l) = acc_p
            Q (n,l) = acc_q
            Qa(n,l) = acc_a
            Qo(n,l) = acc_o
          end do
        end do

!       ---- assemble the components that survive compression
        counter = 0
        do n = 1, n_max
          rn  = radial_exp_coeff(n,k2)
          rdn = radial_exp_coeff_der(n,k2)
          do np = n, n_max
            rnp  = radial_exp_coeff(np,k2)
            rdnp = radial_exp_coeff_der(np,k2)
            base = sq2
            if( n == np ) base = 1.d0
            pre = fourpi*base
            do l = 0, l_max
              counter = counter + 1
              if( skip(counter) )cycle
              this_rad(counter) = pre*( rdn*P(np,l) + rn*Q(np,l) + rdnp*P(n,l) + rnp*Q(n,l) )
              this_azi(counter) = pre*( rn*Qa(np,l) + rnp*Qa(n,l) )
              this_pol(counter) = pre*( rn*Qo(np,l) + rnp*Qo(n,l) )
            end do
          end do
        end do

        if( compress_soap )then
          do k = 1, n_soap
            soap_rad_der(k, k2) = 0.d0
            soap_azi_der(k, k2) = 0.d0
            soap_pol_der(k, k2) = 0.d0
          end do
          do k = 1, compress_P_nonzero
            soap_rad_der(compress_P_i(k), k2) = soap_rad_der(compress_P_i(k), k2) + compress_P_el(k) * &
                                                this_rad(compress_P_j(k))
            soap_azi_der(compress_P_i(k), k2) = soap_azi_der(compress_P_i(k), k2) + compress_P_el(k) * &
                                                this_azi(compress_P_j(k))
            soap_pol_der(compress_P_i(k), k2) = soap_pol_der(compress_P_i(k), k2) + compress_P_el(k) * &
                                                this_pol(compress_P_j(k))
          end do
        else
          soap_rad_der(1:n_soap, k2) = this_rad(1:n_soap)
          soap_azi_der(1:n_soap, k2) = this_azi(1:n_soap)
          soap_pol_der(1:n_soap, k2) = this_pol(1:n_soap)
        end if
      end do
    end do

    deallocate( this_rad, this_azi, this_pol )

  return
  end subroutine soap_derivatives_from_expansion
!**************************************************************************

!**************************************************************************











!**************************************************************************







!**************************************************************************
  subroutine assign_species_multiplicity(max_species_multiplicity, species_types_soap, xyz_species, &
                                         xyz_species_supercell, n_species_soap, all_atoms, which_atom, &
                                         indices, n_neigh, n_atom_pairs, neighbors_list, mask_species, &
                                         species, species_multiplicity, species_multiplicity_supercell )
!                                         species, species_multiplicity )

    implicit none

!   Input variables
    integer, intent(in) :: max_species_multiplicity, n_species_soap, which_atom, indices(1:3), &
                           n_neigh(:), n_atom_pairs, neighbors_list(:)
    character*8, intent(in) :: species_types_soap(:), xyz_species(:), xyz_species_supercell(:)
    logical, intent(in) :: all_atoms

!   Output variables
    integer, allocatable, intent(out) :: species(:,:), species_multiplicity(:), &
                                         species_multiplicity_supercell(:)
!    integer, allocatable :: species_multiplicity_supercell(:)
    logical, allocatable, intent(out) :: mask_species(:,:)

!   Internal variables
    integer, allocatable :: species_supercell(:,:)
    integer :: n_sites, n_sites_supercell, i, j, k, i2, j2, k2, ijunk, counter, i3


    n_sites = size(xyz_species)
    n_sites_supercell = size(xyz_species_supercell)

    allocate( species(1:max_species_multiplicity, 1:n_sites) )
    allocate( species_multiplicity(1:n_sites) )
    species = 0
    species_multiplicity = 0

    do i = 1, n_sites
      do j = 1, n_species_soap
        if( xyz_species(i) == species_types_soap(j) )then
          species_multiplicity(i) = species_multiplicity(i) + 1
          species(species_multiplicity(i), i) = j
        end if
      end do
    end do

    if( n_sites_supercell > n_sites )then
      allocate( species_supercell(1:max_species_multiplicity, 1:n_sites_supercell) )
      allocate( species_multiplicity_supercell(1:n_sites_supercell) )
      species_supercell = 0
      species_multiplicity_supercell = 0
!      counter = 0
!      do i2 = 1, indices(1)
!        do j2 = 1, indices(2)
!          do k2 = 1, indices(3)
!            do i = 1, n_sites
!              counter = counter + 1
!              species_supercell(:, counter) = species(:, i)
!            end do
!          end do
!        end do
!      end do
      do i = 1, n_sites_supercell
        do j = 1, n_species_soap
          if( xyz_species_supercell(i) == species_types_soap(j) )then
            species_multiplicity_supercell(i) = species_multiplicity_supercell(i) + 1
            species_supercell(species_multiplicity_supercell(i), i) = j
          end if
        end do
      end do
    else
      allocate( species_supercell(1:max_species_multiplicity, 1:n_sites_supercell) )
      allocate( species_multiplicity_supercell(1:n_sites_supercell) )
      species_supercell = species
      species_multiplicity_supercell = species_multiplicity
    end if


!   This is perhaps not the most efficient way to select only one atom, fix in the future <----- FIX THIS
    if( .not. all_atoms )then
      n_sites = 1
      deallocate( species )
      allocate( species(1:max_species_multiplicity, 1:n_sites) )
      species(1:max_species_multiplicity, 1) = species_supercell(1:max_species_multiplicity, which_atom)
      species_supercell(1:max_species_multiplicity, which_atom) = species_supercell(1:max_species_multiplicity, 1)
      species_supercell(1:max_species_multiplicity, 1) = species(1:max_species_multiplicity, 1)

      ijunk = species_multiplicity(which_atom)
      deallocate( species_multiplicity )
      allocate( species_multiplicity(1:n_sites) )
      species_multiplicity(1) = ijunk
    end if

    allocate( mask_species(1:n_atom_pairs, 1:n_species_soap) )
    mask_species = .false.

    k2 = 0
    do i = 1, n_sites
      do k = 1, n_neigh(i)
        k2 = k2 + 1
        j = neighbors_list(k2)
        if( k == 1 )then
          i3 = j 
        end if
!        if( k == 1 )then
!          do i2 = 1, species_multiplicity_supercell(i)
!            mask_species(k2, species_supercell(i2, j)) = .true.
!          end do
!        else
!          do i2 = 1, species_multiplicity_supercell(i)
!            mask_species(k2, species_supercell(i2, j)) = .true.
!          end do
!        end if
          do i2 = 1, species_multiplicity_supercell(i3)
            if( species_supercell(i2, j) > 0 )then
              mask_species(k2, species_supercell(i2, j)) = .true.
            end if
          end do
      end do
    end do


!    deallocate( species_supercell, species_multiplicity_supercell )
    deallocate( species_supercell )


  end subroutine
!**************************************************************************


end module soap_turbo_desc
