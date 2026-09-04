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

! ---------------------------------------------------------------------------
! Central-atom second derivatives.
!
! Set soap_hessian_enabled before calling get_soap and it will additionally ask
! the radial expansion for d^2/drj^2 and keep the state that
! get_soap_central_hessian needs. That routine has to be called on the same
! batch, before the next get_soap overwrites the buffers; hs_ready plus the
! recorded batch dimensions are what catches a caller that does not.
!
! The flag costs one extra (n_max, n_atom_pairs) real array and the radial
! kernel's second-derivative branch. It is off by default, so a normal
! energy/force run is untouched.
! ---------------------------------------------------------------------------
  logical, public, save :: soap_hessian_enabled = .false.

! Internal consistency check, off by default. The Hessian rebuilds the angular
! coefficient and its Cartesian gradient from scratch (the Q_l^m route, see
! get_soap_central_hessian), so the same contraction that produces the second
! derivative also produces the first one -- and the first one is already known,
! independently, from soap_rad_der/azi/pol. Comparing them tests the whole
! reconstruction (Omega, dOmega, F, F', the S and Z build, compression and the
! species mask) without any finite differences. It says nothing about the
! second-order pieces, which is exactly what makes it useful for localising a
! failure. The largest relative discrepancy is left in the _err variable.
  logical,  public, save :: soap_hessian_selftest = .false.
  real(dp), public, save :: soap_hessian_selftest_err = 0.d0

  real(dp), allocatable, save, target, private :: wr_rad_der2(:)
  real(dp), allocatable, save, private :: hs_sqrt_dot_p(:)
  logical,  allocatable, save, private :: hs_skip(:)
  logical,  save, private :: hs_ready = .false.
  integer,  save, private :: hs_n_sites = 0, hs_n_pairs = 0, hs_n_max = 0, hs_k_max = 0
  integer,  save, private :: hs_l_max = -1, hs_n_soap = 0

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
  real(dp), pointer, contiguous :: radial_exp_coeff_der2(:,:)
  logical :: do_hessian
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
  character(len=32) :: dump_env
  integer :: dump_len
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

! The second derivatives are only ever wanted alongside the first ones, and the
! radial kernel can only produce them with the corrected filter seed.
  do_hessian = soap_hessian_enabled .and. do_derivatives
  hs_ready = .false.



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
    if( do_hessian )then
      call ensure_real_buffer( wr_rad_der2, int(n_max,i64)*nb_pair )
      radial_exp_coeff_der2(1:n_max, 1:n_atom_pairs)          => wr_rad_der2(1:int(n_max,i64)*nb_pair)
    else
      call ensure_real_buffer( wr_rad_der2, int(n_max,i64) )
      radial_exp_coeff_der2(1:n_max, 1:1)                     => wr_rad_der2(1:int(n_max,i64))
    end if
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
    call ensure_real_buffer( wr_rad_der2, int(n_max,i64) )
    radial_exp_coeff_der2(1:n_max, 1:1) => wr_rad_der2(1:int(n_max,i64))
  end if


  if( do_timing )then
    call cpu_time(time2)
    memory_time = memory_time + time2 - time1
    time1 = time2
  end if


  if( do_hessian .and. basis == "poly3operator" )then
    write(*,*) "ERROR: soap_turbo: the central-atom Hessian needs second radial derivatives, which the"
    write(*,*) "       poly3operator basis does not provide. Use poly3 or poly3gauss."
    stop
  end if

  do i = 1, n_species
    if( basis == "poly3gauss" )then
      if( do_hessian )then
        call get_radial_expansion_coefficients_poly3gauss(n_sites, n_neigh, rjs, alpha_max(i), rcut_soft(i), &
                                                          rcut_hard(i), atom_sigma_r(i), atom_sigma_r_scaling(i), &
                                                          amplitude_scaling(i), nf(i), W(i_beg(i):i_end(i),i_beg(i):i_end(i)), &
                                                          scaling_mode, mask(:,i), radial_enhancement, do_derivatives, &
                                                          radial_exp_coeff(i_beg(i):i_end(i), :), &
                                                          radial_exp_coeff_der(i_beg(i):i_end(i), :), &
                                                          radial_exp_coeff_der2(i_beg(i):i_end(i), :) )
        radial_exp_coeff_der2(i_beg(i):i_end(i), :) = radial_exp_coeff_der2(i_beg(i):i_end(i), :) * global_scaling(i)
      else
        call get_radial_expansion_coefficients_poly3gauss(n_sites, n_neigh, rjs, alpha_max(i), rcut_soft(i), &
                                                          rcut_hard(i), atom_sigma_r(i), atom_sigma_r_scaling(i), &
                                                          amplitude_scaling(i), nf(i), W(i_beg(i):i_end(i),i_beg(i):i_end(i)), &
                                                          scaling_mode, mask(:,i), radial_enhancement, do_derivatives, &
                                                          radial_exp_coeff(i_beg(i):i_end(i), :), &
                                                          radial_exp_coeff_der(i_beg(i):i_end(i), :) )
      end if
      radial_exp_coeff(i_beg(i):i_end(i), :) = radial_exp_coeff(i_beg(i):i_end(i), :) * global_scaling(i)
      radial_exp_coeff_der(i_beg(i):i_end(i), :) = radial_exp_coeff_der(i_beg(i):i_end(i), :) * global_scaling(i)
    else if( basis == "poly3" )then
      if( do_hessian )then
        call get_radial_expansion_coefficients_poly3(n_sites, n_neigh, rjs, alpha_max(i), rcut_soft(i), &
                                                     rcut_hard(i), atom_sigma_r(i), atom_sigma_r_scaling(i), &
                                                     amplitude_scaling(i), nf(i), W(i_beg(i):i_end(i),i_beg(i):i_end(i)), &
                                                     scaling_mode, mask(:,i), radial_enhancement, do_derivatives, &
                                                     do_central(i), central_weight(i), &
                                                     radial_exp_coeff(i_beg(i):i_end(i), :), &
                                                     radial_exp_coeff_der(i_beg(i):i_end(i), :), &
                                                     radial_exp_coeff_der2(i_beg(i):i_end(i), :) )
        radial_exp_coeff_der2(i_beg(i):i_end(i), :) = radial_exp_coeff_der2(i_beg(i):i_end(i), :) * global_scaling(i)
      else
        call get_radial_expansion_coefficients_poly3(n_sites, n_neigh, rjs, alpha_max(i), rcut_soft(i), &
                                                     rcut_hard(i), atom_sigma_r(i), atom_sigma_r_scaling(i), &
                                                     amplitude_scaling(i), nf(i), W(i_beg(i):i_end(i),i_beg(i):i_end(i)), &
                                                     scaling_mode, mask(:,i), radial_enhancement, do_derivatives, &
                                                     do_central(i), central_weight(i), &
                                                     radial_exp_coeff(i_beg(i):i_end(i), :), &
                                                     radial_exp_coeff_der(i_beg(i):i_end(i), :) )
      end if
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

! Debug: dump the radial expansion coefficients and their derivatives at full
! precision, so a GPU build can be compared against this one line for line.
! Enabled by setting TURBOGAP_DUMP_RADIAL; off by default.
  call get_environment_variable( "TURBOGAP_DUMP_RADIAL", dump_env, dump_len )
  if( dump_len > 0 )then
    open(unit=78, file="radial_exp_coeff_dump.dat", status="unknown", position="append")
    write(78,'(A,2I8)') "# n_atom_pairs n_max ", n_atom_pairs, n_max
!   With do_derivatives off, radial_exp_coeff_der is a one-column dummy -- see
!   its allocation above -- so it is only safe to index for k2 > 1 when the
!   derivatives were actually computed. Write a zero in that column instead of
!   reading off the end of the array.
    do k2 = 1, n_atom_pairs
      do i = 1, n_max
        if( do_derivatives )then
          write(78,'(2I8,1X,ES26.16,1X,ES26.16,1X,ES26.16)') k2, i, rjs(k2), &
            radial_exp_coeff(i,k2), radial_exp_coeff_der(i,k2)
        else
          write(78,'(2I8,1X,ES26.16,1X,ES26.16,1X,ES26.16)') k2, i, rjs(k2), &
            radial_exp_coeff(i,k2), 0.d0
        end if
      end do
    end do
    close(78)
  end if



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
! Everything get_soap_central_hessian needs that is not either a persistent
! buffer or an argument the caller already holds. sqrt_dot_p is the one piece
! of the normalisation that cannot be recovered from the returned soap, and the
! skip mask decides which uncompressed components a weight vector may touch.
  if( do_hessian )then
    if( allocated(hs_sqrt_dot_p) )then
      if( size(hs_sqrt_dot_p) < n_sites ) deallocate( hs_sqrt_dot_p )
    end if
    if( .not. allocated(hs_sqrt_dot_p) ) allocate( hs_sqrt_dot_p(1:n_sites) )
    hs_sqrt_dot_p(1:n_sites) = sqrt_dot_p(1:n_sites)
    if( allocated(hs_skip) )then
      if( size(hs_skip) /= size(skip_soap_component_flattened) ) deallocate( hs_skip )
    end if
    if( .not. allocated(hs_skip) ) allocate( hs_skip(1:size(skip_soap_component_flattened)) )
    hs_skip = skip_soap_component_flattened
    hs_n_sites = n_sites
    hs_n_pairs = n_atom_pairs
    hs_n_max   = n_max
    hs_k_max   = k_max
    hs_l_max   = l_max
    hs_n_soap  = n_soap
    hs_ready   = .true.
  end if

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




!**************************************************************************
!
! d^2 soap / d(central atom) d(neighbour), contracted against weight vectors.
!
! WHAT THIS IS FOR
!
! A GAP dipole model predicts mu_i = dE_i/dr_i, the local energy differentiated
! with respect to its own central atom. Its Cartesian gradient is therefore a
! second derivative of the descriptor, and the piece that soap_turbo has to
! supply is
!
!     d^2 q_d / d r_ia d r_jb
!
! for site i, one of its neighbours j, and every descriptor component d. The
! full Hessian -- every (j,k) pair of neighbours -- is O(n_neigh^2) blocks per
! site and is not wanted here: only the row in which one index is the central
! atom. That row is O(n_neigh), and getting it for that price is the whole
! point of this routine.
!
! HOW THE CENTRAL-ATOM ROW COLLAPSES
!
! The descriptor depends on positions only through the displacements
! u_j = r_j - r_i, so d/dr_ia = -sum_k d/du_ka and
!
!     d^2 q / d r_ia d r_jb  =  - sum_k  d^2 q / d u_ka d u_jb.            (1)
!
! The power spectrum is quadratic in the expansion coefficients,
! p_(n,np,l) = sum_m mult * Re( c_n conjg(c_np) ), and c is a plain sum over
! neighbours, so d^2 c / du_k du_j vanishes unless k = j. Writing
!
!     D^j_(n,k,b) = d c_(n,k) / d u_jb            (one pair, first derivative)
!     E^j_(n,k,ab) = d^2 c_(n,k) / d u_ja d u_jb  (one pair, 3x3 Hessian)
!     G_(n,k,a)   = sum_j D^j_(n,k,a)             (accumulated over the site)
!
! the sum over k in (1) telescopes into G, and only the j = k term keeps a
! second derivative:
!
!     d^2 p / dr_ia dr_jb = - sum_m mult Re[ E^j_(n,ab) conjg(c_np)
!                                          + c_n conjg(E^j_(np,ab))
!                                          + G_(n,a) conjg(D^j_(np,b))
!                                          + D^j_(n,b) conjg(G_(np,a)) ].
!
! So one accumulator per site replaces the sum over every other neighbour, and
! the per-pair cost stays proportional to n_max*k_max. G is the same object the
! first-derivative code already forms as the self pair's soap_cart_der.
!
! CONTRACTING EARLY
!
! What a dipole gradient needs is not the Hessian of each descriptor component
! but sum_d w_d d^2 q_d / dr_ia dr_jb with w = dE_i/dq. Storing the
! per-component Hessian would be (3,3,n_soap,n_atom_pairs) -- 144 GB for a
! 7000-atom cell -- so the weight vector is pushed all the way inside instead
! and the output is (3,3,n_vec,n_atom_pairs). With the symmetric matrix
!
!     S_(n,np,l) = 2 v_(n,n,l)          n = np
!                = sqrt(2) v_(n,np,l)   n /= np      (both orders)
!
! built once per site from the weight vector, and Z_(n,k) = sum_np S c_(np,k),
! ZG_(n,k,a) = sum_np S G_(np,k,a), the two pairs of terms above collapse to
!
!     v . d^2 p / dr_ia dr_jb = - sum_k w(m) { sum_n Re[ E^j_(n,ab) conjg(Z_n) ]
!                                            + sum_n Re[ D^j_(n,b) conjg(ZG_(n,a)) ] }
!
! Both pairs of terms carry no extra factor: S already folds in the 2 on the
! diagonal and the sqrt(2) off it, which is precisely what turns a sum over
! n <= np of "X_n conjg(c_np) + c_n conjg(X_np)" into a full double sum over
! one ordering. Compression is a linear map, so a compressed weight vector is
! handled by transforming it once with P^T before S is built.
!
! THE SINGLE-PAIR HESSIAN, WITHOUT THE POLES
!
! E^j still needs the second Cartesian derivative of one pair's contribution,
! c = g_(nl)(r) * Omega_(lm)(nhat). Doing that through r, theta, phi produces
! 1/sin(theta) and 1/sin^2(theta) factors that cancel only in the total, which
! is how the first-derivative code ended up returning Plm/sin and sin*dPlm/dx.
! There is no comfortable second-derivative analogue of that trick.
!
! Instead the angular factor is written in a form that has no polar singularity
! at all. With the Condon-Shortley convention used by get_plm_array,
! P_l^m(x) = (-1)^m (1-x^2)^(m/2) Q_l^m(x) where Q_l^m = d^m P_l / dx^m is a
! plain polynomial, and with t = z/r and w = (x - iy)/r,
!
!     P_l^m(cos t) e^(-i m phi) = (-1)^m Q_l^m(t) w^m,
!
! so Omega = C_k Q_l^m(t) w^m is a smooth function of the Cartesian
! displacement everywhere except r = 0. Its gradient and Hessian follow from
! those of t and w, which are elementary. Q, Q' and Q'' come from one table:
! Q' = Q_l^(m+1) and Q'' = Q_l^(m+2), zero above m = l.
!
! COST
!
! Per pair: n_max*k_max complex multiply-adds for G, and 9*n_max*k_max for the
! contraction, per weight vector. Roughly three times the existing force
! contraction, and linear in the number of pairs.
!
! USAGE
!
! Set soap_hessian_enabled = .true., call get_soap with do_derivatives, then
! call this on the same batch before the next get_soap overwrites the buffers.
! It needs legacy_filter_seed = .false. (the radial second derivatives do) and
! a poly3 or poly3gauss basis.
!
! hess(:,:,v,k2) for pair k2 = (site i, neighbour j) is
! sum_d vecs(d,v,i) d^2 q_d/dr_ia dr_jb, with the first index a and the second
! b. It is not symmetric: a belongs to the central atom and b to the
! neighbour. For j = 1, the self pair, the block returned is the central
! atom's own, d^2/dr_ia dr_ib, built as minus the sum of the others.
!
  subroutine get_soap_central_hessian(n_sites, n_neigh, n_atom_pairs, n_species, mask, &
                                      rjs, thetas, phis, l_max, rcut_hard, &
                                      atom_sigma_t, atom_sigma_t_scaling, &
                                      compress_soap, compress_P_nonzero, compress_P_i, &
                                      compress_P_j, compress_P_el, soap, n_vec, vecs, hess)

    implicit none

!   Input
    integer,  intent(in) :: n_sites, n_neigh(:), n_atom_pairs, n_species, l_max, n_vec
    integer,  intent(in) :: compress_P_nonzero, compress_P_i(:), compress_P_j(:)
    logical,  intent(in) :: mask(:,:), compress_soap
    real(dp), intent(in) :: rjs(:), thetas(:), phis(:), rcut_hard(:)
    real(dp), intent(in) :: atom_sigma_t(:), atom_sigma_t_scaling(:), compress_P_el(:)
!   soap as get_soap returned it, i.e. normalised
    real(dp), intent(in) :: soap(:,:)
!   vecs(1:n_soap, 1:n_vec, 1:n_sites): the weight vectors, in the same
!   (normalised, compressed) space as soap
    real(dp), intent(in) :: vecs(:,:,:)

!   Output
    real(dp), intent(out) :: hess(:,:,:,:)

!   Views of the buffers get_soap filled
    complex(dp), pointer, contiguous :: cnk(:,:,:)
    real(dp),    pointer, contiguous :: rad(:,:), rad_d(:,:), rad_d2(:,:)
    real(dp),    pointer, contiguous :: srad(:,:), sazi(:,:), spol(:,:)

!   Per-site
    complex(dp), allocatable :: Gacc(:,:,:), Zv(:,:,:), ZGv(:,:,:,:)
    real(dp),    allocatable :: Smat(:,:,:), vunc(:), pvec(:), dp1(:,:), dp2(:,:)
    real(dp),    allocatable :: hv(:,:,:), hsum(:,:,:), u2rec(:,:)
    real(dp) :: refd, chk_err, chk_ref

!   Per-pair angular factors
    complex(dp), allocatable :: Om(:), dOm(:,:), d2Om(:,:,:)
    complex(dp) :: wc, wvec(3), whes(3,3), ec(3)
    complex(dp) :: zc, gc, a0, a1, a2, b0(3), b1(3), z0, z1, z2, gg0(3), gg1(3)
    complex(dp) :: tE, tD, om_k, dom_k(3)
    complex(dp), allocatable :: wpow(:)

!   Tables
    real(dp), allocatable :: preflm(:), Ck(:), Qt(:), wm(:), fact(:)
    real(dp), allocatable :: il(:), dilz(:), d2ilz(:), Fl(:), Fl1(:), Fl2(:)
    integer,  allocatable :: lofk(:), mofk(:)

    real(dp) :: nh(3), tvec(3), thes(3,3), pi, fourpi, rcut_max, rr, invr, invr2
    real(dp) :: sig, sscal, amp, ampd, ampdd, zz, invz, zprime, zprime2, xarg
    real(dp) :: sin_t, cos_t, sin_f, cos_f, cxr, cyr, czr, pxr, pyr, pzr, axr, ayr
    real(dp) :: tt, f0, f1, f2, r0, r1, r2, sqrt2, wmk, grad0, grad1
    real(dp) :: s, invs, invs3, invs5, wp, u1(3), u2(3), t1(3), t2(3), c12(3,3)
    real(dp) :: q0, q1, q2

    integer :: n_soap, n_soap_unc, n_max, k_max, nv, i, j, k, k2, k3, n, np, l, m
    integer :: a, b, v, counter, kb, i_sp, isp_found
    integer(i64) :: nb_pair
    logical :: ok

!   ---------------------------------------------------------------------
!   State check. Everything below reads buffers that the next get_soap call
!   will overwrite, so a caller that gets the order wrong has to be told.
!   ---------------------------------------------------------------------
    if( .not. hs_ready )then
      write(*,*) "ERROR: get_soap_central_hessian was called without a matching get_soap."
      write(*,*) "       Set soap_hessian_enabled = .true. and call get_soap with"
      write(*,*) "       do_derivatives = .true. first."
      stop
    end if
    n_soap = size(soap, 1)
    if( n_sites /= hs_n_sites .or. n_atom_pairs /= hs_n_pairs .or. &
        l_max /= hs_l_max .or. n_soap /= hs_n_soap )then
      write(*,*) "ERROR: get_soap_central_hessian was called on a different batch than the"
      write(*,*) "       get_soap that filled the buffers it reads."
      stop
    end if

    pi     = dacos(-1.d0)
    fourpi = 4.d0*pi
    sqrt2  = dsqrt(2.d0)
    n_max  = hs_n_max
    k_max  = hs_k_max
    n_soap_unc = size(hs_skip)
    rcut_max = maxval(rcut_hard)
    nv = n_vec + 1                       ! the user's vectors, plus p itself

    nb_pair = int(n_atom_pairs, i64)
    cnk   (1:k_max, 1:n_max, 1:n_sites) => wc_cnk     (1:int(k_max,i64)*int(n_max,i64)*int(n_sites,i64))
    rad   (1:n_max, 1:n_atom_pairs)     => wr_rad     (1:int(n_max,i64)*nb_pair)
    rad_d (1:n_max, 1:n_atom_pairs)     => wr_rad_der (1:int(n_max,i64)*nb_pair)
    rad_d2(1:n_max, 1:n_atom_pairs)     => wr_rad_der2(1:int(n_max,i64)*nb_pair)
    srad  (1:n_soap, 1:n_atom_pairs)    => wr_soap_rad(1:int(n_soap,i64)*nb_pair)
    sazi  (1:n_soap, 1:n_atom_pairs)    => wr_soap_azi(1:int(n_soap,i64)*nb_pair)
    spol  (1:n_soap, 1:n_atom_pairs)    => wr_soap_pol(1:int(n_soap,i64)*nb_pair)

!   ---------------------------------------------------------------------
!   (l,m) tables. k = 1 + l*(l+1)/2 + m, so l and m are recovered once here
!   rather than carried through nested loops.
!   ---------------------------------------------------------------------
    allocate( lofk(1:k_max), mofk(1:k_max), wm(1:k_max), preflm(1:k_max), Ck(1:k_max), Qt(1:k_max) )
    k = 0
    do l = 0, l_max
      do m = 0, l
        k = k + 1
        lofk(k) = l
        mofk(k) = m
        if( m > 0 )then
          wm(k) = 2.d0
        else
          wm(k) = 1.d0
        end if
      end do
    end do
    call get_preflm(preflm, l_max)
!   Omega = Ck * Q_l^m(t) * w^m carries the (-1)^m that turns P_l^m into Q_l^m
    do k = 1, k_max
      if( mod(mofk(k), 2) == 0 )then
        Ck(k) =  preflm(k)
      else
        Ck(k) = -preflm(k)
      end if
    end do

    allocate( fact(1:max(l_max+2,1)) )
    call get_semifactorial_array(fact, l_max+2)
    allocate( il(0:l_max+2), dilz(0:l_max), d2ilz(0:l_max) )
    allocate( Fl(0:l_max), Fl1(0:l_max), Fl2(0:l_max) )
    allocate( Om(1:k_max), dOm(1:3,1:k_max), d2Om(1:3,1:3,1:k_max) )
    allocate( wpow(-2:max(l_max,2)) )

    allocate( Gacc(1:n_max, 1:k_max, 1:3) )
    allocate( Zv  (1:n_max, 1:k_max, 1:nv) )
    allocate( ZGv (1:n_max, 1:k_max, 1:3, 1:nv) )
    allocate( Smat(1:n_max, 1:n_max, 0:l_max) )
    allocate( vunc(1:n_soap_unc), pvec(1:n_soap) )
    allocate( dp1(1:n_soap, 1:3), dp2(1:n_soap, 1:3) )
    allocate( hv(1:3, 1:3, 1:nv), hsum(1:3, 1:3, 1:n_vec) )
    allocate( u2rec(1:3, 1:nv) )

    hess = 0.d0
    chk_err = 0.d0
    chk_ref = 0.d0

!   =====================================================================
    k3 = 0
    k2 = 0
    do i = 1, n_sites
      s = hs_sqrt_dot_p(i)
      invs  = 1.d0/s
      invs3 = invs/s**2
      invs5 = invs3/s**2
      pvec(1:n_soap) = soap(1:n_soap, i) * s

!     ---- pass A: G and the unnormalised central-atom derivative ------
      Gacc = (0.d0, 0.d0)
      dp1  = 0.d0
      k3   = k2 + 1
      do j = 1, n_neigh(i)
        k2 = k2 + 1
        if( j == 1 ) cycle
        call pair_geometry(k2, ok, .false.)
        if( .not. ok ) cycle
        call pair_dp(k2, dp2)
        do a = 1, 3
          do k = 1, n_soap
            dp1(k, a) = dp1(k, a) - dp2(k, a)
          end do
        end do
        do k = 1, k_max
          l = lofk(k)
          f0 = Fl(l)
          f1 = Fl1(l)
          om_k = Om(k)
          dom_k = dOm(1:3, k)
          do n = 1, n_max
            r0 = fourpi*rad(n, k2)
            r1 = fourpi*rad_d(n, k2)
            grad0 = r0*f0                    ! g
            grad1 = r1*f0 + r0*f1            ! g'
            do a = 1, 3
              Gacc(n, k, a) = Gacc(n, k, a) + grad1*nh(a)*om_k + grad0*dom_k(a)
            end do
          end do
        end do
      end do

!     ---- the weight vectors, uncompressed, and their S matrices ------
      do v = 1, nv
        if( v <= n_vec )then
          call spread_weights(vecs(1:n_soap, v, i), vunc)
        else
          call spread_weights(pvec, vunc)
        end if
        Smat = 0.d0
        counter = 0
        do n = 1, n_max
          do np = n, n_max
            do l = 0, l_max
              counter = counter + 1
              if( n == np )then
                Smat(n, n, l) = 2.d0*vunc(counter)
              else
                Smat(n, np, l) = sqrt2*vunc(counter)
                Smat(np, n, l) = Smat(n, np, l)
              end if
            end do
          end do
        end do
        do k = 1, k_max
          l = lofk(k)
          do n = 1, n_max
            zc = (0.d0, 0.d0)
            b0 = (0.d0, 0.d0)
            do np = 1, n_max
              zc = zc + Smat(n, np, l)*cnk(k, np, i)
              do a = 1, 3
                b0(a) = b0(a) + Smat(n, np, l)*Gacc(np, k, a)
              end do
            end do
            Zv(n, k, v) = zc
            do a = 1, 3
              ZGv(n, k, a, v) = b0(a)
            end do
          end do
        end do
      end do

!     ---- pass B: the per-pair contraction ---------------------------
      hsum = 0.d0
      k2 = k3 - 1
      do j = 1, n_neigh(i)
        k2 = k2 + 1
        if( j == 1 ) cycle
        call pair_geometry(k2, ok, .true.)
        if( .not. ok ) cycle
        call pair_dp(k2, dp2)

        hv = 0.d0
        u2rec = 0.d0
        do k = 1, k_max
          l   = lofk(k)
          f0  = Fl(l)
          f1  = Fl1(l)
          f2  = Fl2(l)
          wmk = wm(k)
          om_k  = Om(k)
          dom_k = dOm(1:3, k)
          do v = 1, nv
            a0 = (0.d0, 0.d0); a1 = (0.d0, 0.d0); a2 = (0.d0, 0.d0)
            b0 = (0.d0, 0.d0); b1 = (0.d0, 0.d0)
            do n = 1, n_max
              r0 = rad(n, k2); r1 = rad_d(n, k2); r2 = rad_d2(n, k2)
              zc = dconjg(Zv(n, k, v))
              a0 = a0 + r0*zc
              a1 = a1 + r1*zc
              a2 = a2 + r2*zc
              do a = 1, 3
                gc = dconjg(ZGv(n, k, a, v))
                b0(a) = b0(a) + r0*gc
                b1(a) = b1(a) + r1*gc
              end do
            end do
            z0 = fourpi*f0*a0
            z1 = fourpi*(f0*a1 + f1*a0)
            z2 = fourpi*(f0*a2 + 2.d0*f1*a1 + f2*a0)
            do a = 1, 3
              gg0(a) = fourpi*f0*b0(a)
              gg1(a) = fourpi*(f0*b1(a) + f1*b0(a))
            end do
!           The first derivative falls out of the same z0/z1: it is
!           v . dp/du_jb = sum_k w(m) Re[ z1 n_b Omega + z0 dOmega_b ].
            if( soap_hessian_selftest )then
              do b = 1, 3
                u2rec(b, v) = u2rec(b, v) + wmk*dble( z1*nh(b)*om_k + z0*dom_k(b) )
              end do
            end if
            do b = 1, 3
              do a = 1, 3
                tE = (z2 - z1*invr)*nh(a)*nh(b)*om_k &
                     + z1*(nh(a)*dom_k(b) + nh(b)*dom_k(a)) &
                     + z0*d2Om(a, b, k)
                if( a == b ) tE = tE + z1*invr*om_k
                tD = gg1(a)*nh(b)*om_k + gg0(a)*dom_k(b)
                hv(a, b, v) = hv(a, b, v) - wmk*( dble(tE) + dble(tD) )
              end do
            end do
          end do
        end do

        if( soap_hessian_selftest )then
          do v = 1, nv
            do b = 1, 3
              if( v <= n_vec )then
                refd = dot_product(vecs(1:n_soap, v, i), dp2(:, b))
              else
                refd = dot_product(pvec, dp2(:, b))
              end if
              chk_err = max( chk_err, abs( u2rec(b, v) - refd ) )
              chk_ref = max( chk_ref, abs( refd ) )
            end do
          end do
        end if

!       ---- normalisation ------------------------------------------
!       q = p/s, so with d1 = d/dr_ia and d2 = d/dr_jb,
!         w.d1d2 q = (w.d1d2 p)/s - [(w.d1 p)(p.d2 p) + (w.d2 p)(p.d1 p)]/s^3
!                    - (w.p)[(d1 p . d2 p) + (p.d1d2 p)]/s^3
!                    + 3 (w.p)(p.d1 p)(p.d2 p)/s^5
!       The p-weighted Hessian it needs is hv(:,:,nv), which is why p was
!       carried through the contraction as an extra weight vector.
        do a = 1, 3
          t1(a) = dot_product(pvec, dp1(:, a))
          t2(a) = dot_product(pvec, dp2(:, a))
        end do
        do b = 1, 3
          do a = 1, 3
            c12(a, b) = dot_product(dp1(:, a), dp2(:, b))
          end do
        end do
        do v = 1, n_vec
          wp = dot_product(vecs(1:n_soap, v, i), pvec)
          do a = 1, 3
            u1(a) = dot_product(vecs(1:n_soap, v, i), dp1(:, a))
            u2(a) = dot_product(vecs(1:n_soap, v, i), dp2(:, a))
          end do
          do b = 1, 3
            do a = 1, 3
              hess(a, b, v, k2) = hv(a, b, v)*invs &
                                  - ( u1(a)*t2(b) + u2(b)*t1(a) )*invs3 &
                                  - wp*( c12(a, b) + hv(a, b, nv) )*invs3 &
                                  + 3.d0*wp*t1(a)*t2(b)*invs5
              hsum(a, b, v) = hsum(a, b, v) + hess(a, b, v, k2)
            end do
          end do
        end do
      end do

!     The central atom's own block. A rigid shift of every atom leaves the
!     descriptor alone, so the blocks of one row sum to zero.
      do v = 1, n_vec
        do b = 1, 3
          do a = 1, 3
            hess(a, b, v, k3) = -hsum(a, b, v)
          end do
        end do
      end do
    end do

    if( soap_hessian_selftest )then
      soap_hessian_selftest_err = chk_err / max( chk_ref, 1.d-300 )
    end if

    deallocate( lofk, mofk, wm, preflm, Ck, Qt, fact, il, dilz, d2ilz, Fl, Fl1, Fl2 )
    deallocate( Om, dOm, d2Om, wpow, Gacc, Zv, ZGv, Smat, vunc, pvec, dp1, dp2, hv, hsum, u2rec )

  contains

!   -------------------------------------------------------------------
!   P^T applied to a weight vector, so that a compressed vector is turned
!   into the uncompressed (n,np,l) space the S matrix is built in. Skipped
!   components are cleared: they are not part of the descriptor and a weight
!   must not reach them.
    subroutine spread_weights(vin, vout)
      real(dp), intent(in)  :: vin(:)
      real(dp), intent(out) :: vout(:)
      integer :: kk
      vout = 0.d0
      if( compress_soap )then
        do kk = 1, compress_P_nonzero
          vout(compress_P_j(kk)) = vout(compress_P_j(kk)) + compress_P_el(kk)*vin(compress_P_i(kk))
        end do
      else
        vout(1:size(vin)) = vin
      end if
      do kk = 1, size(vout)
        if( hs_skip(kk) ) vout(kk) = 0.d0
      end do
    end subroutine spread_weights

!   -------------------------------------------------------------------
!   The unnormalised Cartesian derivative of the compressed SOAP vector for
!   one pair: the same rotation get_soap applies, without the normalisation.
    subroutine pair_dp(kk, d)
      integer, intent(in) :: kk
      real(dp), intent(out) :: d(:,:)
      integer :: kq
      sin_t = dsin(thetas(kk)); cos_t = dcos(thetas(kk))
      sin_f = dsin(phis(kk));   cos_f = dcos(phis(kk))
      invr = 1.d0/rjs(kk)
      cxr =  sin_t*cos_f;        cyr =  sin_t*sin_f;        czr = cos_t
      pxr = -cos_t*cos_f*invr;   pyr = -cos_t*sin_f*invr;   pzr = sin_t*invr
      axr = -sin_f*invr;         ayr =  cos_f*invr
      do kq = 1, n_soap
        d(kq, 1) = cxr*srad(kq, kk) + pxr*spol(kq, kk) + axr*sazi(kq, kk)
        d(kq, 2) = cyr*srad(kq, kk) + pyr*spol(kq, kk) + ayr*sazi(kq, kk)
        d(kq, 3) = czr*srad(kq, kk) + pzr*spol(kq, kk)
      end do
    end subroutine pair_dp

!   -------------------------------------------------------------------
!   Everything about one pair that does not involve the radial basis: the
!   unit vector, the radial part of the angular coefficient (F and its two
!   derivatives) and Omega with its Cartesian gradient and, if asked, its
!   Cartesian Hessian.
    subroutine pair_geometry(kk, good, want2)
      integer, intent(in) :: kk
      logical, intent(out) :: good
      logical, intent(in) :: want2
      integer :: ll, mm, kq, aa, bb, kbb
      real(dp) :: rl

      good = .false.
      rr = rjs(kk)
      if( rr >= rcut_max .or. rr < 1.d-10 ) return
      isp_found = 0
      do i_sp = 1, n_species
        if( mask(kk, i_sp) )then
          isp_found = i_sp
          exit
        end if
      end do
      if( isp_found == 0 ) return
      good = .true.

      invr  = 1.d0/rr
      invr2 = invr*invr

!     ---- F_l(r) = amplitude(r) * ilexp(l, r/sigma) -------------------
      sscal = atom_sigma_t_scaling(isp_found)
      sig   = atom_sigma_t(isp_found) + sscal*rr
      amp   = rcut_max**2 / sig**2
      ampd  = -2.d0*sscal/sig * amp
      ampdd =  6.d0*sscal**2/sig**2 * amp
      xarg  = rr/sig
      zz    = xarg*xarg
      invz  = 1.d0/zz
      call get_ilexp(il, fact, l_max+2, xarg)
!     ilexp(l) = i_l(z) exp(-z) with z = (r/sigma)^2. The upward form of
!     d/dz, D_l = I_(l+1) + (l/z - 1) I_l, keeps every term of the same order
!     in z as I_l itself, so it does not lose digits the way the downward
!     recursion does when z is small.
      do l = 0, l_max
        rl = dfloat(l)
        dilz(l)  = il(l+1) + (rl*invz - 1.d0)*il(l)
        d2ilz(l) = il(l+2) + (dfloat(2*l+1)*invz - 2.d0)*il(l+1) &
                 + (rl*(rl-1.d0)*invz*invz - 2.d0*rl*invz + 1.d0)*il(l)
      end do
      zprime  = 2.d0*rr/sig**2 - 2.d0*rr**2*sscal/sig**3
      zprime2 = 2.d0/sig**2 - 8.d0*rr*sscal/sig**3 + 6.d0*rr**2*sscal**2/sig**4
      do l = 0, l_max
        f0 = il(l)
        f1 = dilz(l)*zprime
        f2 = d2ilz(l)*zprime**2 + dilz(l)*zprime2
        Fl (l) = amp*f0
        Fl1(l) = ampd*f0 + amp*f1
        Fl2(l) = ampdd*f0 + 2.d0*ampd*f1 + amp*f2
      end do

!     ---- geometry ---------------------------------------------------
      sin_t = dsin(thetas(kk)); cos_t = dcos(thetas(kk))
      sin_f = dsin(phis(kk));   cos_f = dcos(phis(kk))
      nh(1) = sin_t*cos_f
      nh(2) = sin_t*sin_f
      nh(3) = cos_t
      tt = nh(3)
      wc = dcmplx(nh(1), -nh(2))
      ec(1) = (1.d0, 0.d0)
      ec(2) = (0.d0, -1.d0)
      ec(3) = (0.d0, 0.d0)

      do aa = 1, 3
        tvec(aa) = -tt*nh(aa)*invr
        wvec(aa) = -wc*nh(aa)*invr
      end do
      tvec(3) = tvec(3) + invr
      do aa = 1, 3
        wvec(aa) = wvec(aa) + ec(aa)*invr
      end do
      if( want2 )then
        do bb = 1, 3
          do aa = 1, 3
            thes(aa, bb) = ( 3.d0*tt*nh(aa)*nh(bb) )*invr2
            whes(aa, bb) = ( 3.d0*wc*nh(aa)*nh(bb) - ec(aa)*nh(bb) - ec(bb)*nh(aa) )*invr2
          end do
        end do
        do aa = 1, 3
          thes(3, aa) = thes(3, aa) - nh(aa)*invr2
          thes(aa, 3) = thes(aa, 3) - nh(aa)*invr2
          thes(aa, aa) = thes(aa, aa) - tt*invr2
          whes(aa, aa) = whes(aa, aa) - wc*invr2
        end do
      end if

!     ---- Q_l^m(t) = d^m P_l/dt^m, and the powers of w ----------------
      Qt(1) = 1.d0
      if( l_max >= 1 )then
        Qt(2) = tt
        Qt(3) = 1.d0
      end if
      do ll = 2, l_max
        kbb = 1 + ll*(ll+1)/2
        do mm = 0, ll-2
          kq = kbb + mm
          Qt(kq) = ( dfloat(2*ll-1)*tt*Qt(kq-ll) - dfloat(ll+mm-1)*Qt(kq-2*ll+1) ) / dfloat(ll-mm)
        end do
        kq = kbb + ll - 1
        Qt(kq) = dfloat(2*ll-1)*tt*Qt(kq-ll)
        kq = kbb + ll
        Qt(kq) = dfloat(2*ll-1)*Qt(kq-ll-1)
      end do
      wpow(-2) = (0.d0, 0.d0)
      wpow(-1) = (0.d0, 0.d0)
      wpow(0)  = (1.d0, 0.d0)
      do mm = 1, l_max
        wpow(mm) = wpow(mm-1)*wc
      end do

!     ---- Omega and its Cartesian derivatives -------------------------
      do kq = 1, k_max
        ll = lofk(kq)
        mm = mofk(kq)
        q0 = Qt(kq)
        if( mm+1 <= ll )then
          q1 = Qt(kq+1)
        else
          q1 = 0.d0
        end if
        if( mm+2 <= ll )then
          q2 = Qt(kq+2)
        else
          q2 = 0.d0
        end if
        Om(kq) = Ck(kq)*q0*wpow(mm)
        do aa = 1, 3
          dOm(aa, kq) = Ck(kq)*( q1*tvec(aa)*wpow(mm) + q0*dfloat(mm)*wpow(mm-1)*wvec(aa) )
        end do
        if( want2 )then
          do bb = 1, 3
            do aa = 1, 3
              d2Om(aa, bb, kq) = Ck(kq)*( q2*tvec(aa)*tvec(bb)*wpow(mm) &
                                        + q1*thes(aa, bb)*wpow(mm) &
                                        + q1*dfloat(mm)*wpow(mm-1)*( tvec(aa)*wvec(bb) + tvec(bb)*wvec(aa) ) &
                                        + q0*dfloat(mm*(mm-1))*wpow(mm-2)*wvec(aa)*wvec(bb) &
                                        + q0*dfloat(mm)*wpow(mm-1)*whes(aa, bb) )
            end do
          end do
        end if
      end do
    end subroutine pair_geometry

  end subroutine get_soap_central_hessian


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
