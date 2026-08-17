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


module soap_turbo_radial

! Kind parameters, declared locally on purpose: soap_turbo is a submodule with
! its own upstream and must not depend on the parent tree's kinds module.
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64


  use soap_turbo_functions

! private so a module that uses another does not inherit a second copy
  integer, parameter, private :: sp  = real32
  integer, parameter, private :: dp  = real64
  integer, parameter, private :: i32 = int32
  integer, parameter, private :: i64 = int64
  private :: real32, real64, int32, int64


! ---------------------------------------------------------------------------
! Filter seed: which of the two filter recursions the radial expansion uses.
!
! The Gaussian filter integral over [rcut_soft, rcut_hard] is built by a
! recursion in alpha that is seeded at alpha = -2. Every step of it drops the
! surface term at r = rcut_hard, on the grounds that (1 - r/rcut_hard)**(alpha+2)
! vanishes there -- true for alpha >= -1, but not for the alpha = -2 seed,
! where the exponent is 0 and the factor is 1.
!
! Keeping that term makes exp_coeff_der agree with a central difference of
! exp_coeff to machine precision (verified in real128). Dropping it leaves a
! relative disagreement of order exp(-nf**2/2) -- 3.4e-4 for the default
! nf = 4 -- concentrated on neighbours sitting at the hard cutoff, which for
! alpha_max = 8 shows up as a ~5e-6 relative error in both the coefficients
! and their derivatives.
!
!   .true.   bit-compatible with the pre-2026 code (default)
!   .false.  finite-difference consistent
!
! Flipping this changes descriptor values at the 1e-5 level near rcut_hard, so
! it changes energies and forces from a fitted GAP by a comparable relative
! amount. Refit or re-baseline before turning it off in production.
! ---------------------------------------------------------------------------
  logical, public :: legacy_filter_seed = .true.




  contains

!**************************************************************************
!
! This function returns the normalization coefficients for the polynomial
! basis functions used to construct the orthonormal basis.
!
  function N_a(rcut, a)
    implicit none

    integer :: a, b
    real(dp) :: rcut, N_a

    b = 2*a + 5

!   **************** New basis ******************
!    N_a = dsqrt( rcut**b / dfloat(b) )
    N_a = dsqrt( rcut / dfloat(b) )
!   *********************************************

  return
  end function
!**************************************************************************









!**************************************************************************
!
! This subroutine returns the radial expansion coefficients using the
! polynomial basis set and the saparable Gaussian representation for the
! atomic sites.
!
! It is a thin wrapper around radial_expansion_kernel(), which is shared
! with the poly3gauss basis; see the comments there.
!
  subroutine get_radial_expansion_coefficients_poly3(n_sites, n_neigh, rjs_in, alpha_max, rcut_soft_in, &
                                                     rcut_hard_in, atom_sigma_in, atom_sigma_scaling, &
                                                     amplitude_scaling, nf, W, scaling_mode, mask, &
                                                     radial_enhancement, do_derivatives, do_central, &
                                                     central_weight, exp_coeff, exp_coeff_der, &
                                                     exp_coeff_der2)
!   Expansion coefficients using the polynomial basis with smooth filter

    implicit none

    integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
    real(dp), intent(in) :: rcut_soft_in, rcut_hard_in, rjs_in(:), atom_sigma_in, nf, atom_sigma_scaling
    real(dp), intent(in) :: amplitude_scaling, central_weight
    logical, intent(in) :: mask(:), do_derivatives, do_central
    character(*), intent(in) :: scaling_mode
    real(dp) :: W(:,:)
!   Results will be stored in exp_coeff, which is an array of dimension
!   (alpha_max, n_atom_pairs)
    real(dp) :: exp_coeff(:,:), exp_coeff_der(:,:)
!   Optional second derivative d^2 exp_coeff / d rj^2, same shape. Requesting it
!   requires do_derivatives and legacy_filter_seed = .false.; see the kernel.
    real(dp), optional :: exp_coeff_der2(:,:)

    logical, save :: print_message = .true.

!   For this basis numerical instabilities in the orthonormal basis construction develop
!   above alpha_max = 7. These are relatively small up to alpha_max = 10 and become
!   catastrophic at alpha_max = 12.
    if( alpha_max > 10 )then
      write(*,*) "-------------------------------------------------------------------------------"
      write(*,*) "ERROR: Due to numerical instabilities in the basis construction for the poly3   <---- ERROR"
      write(*,*) "basis, it is strongly recommended not to exceed alpha_max = 7. For"
      write(*,*) "alpha_max > 10 the instabilities are too large to proceed. You can do your own"
      write(*,*) "testing for 7 < alpha_max < 11, the results might still be useful within that"
      write(*,*) "range. Note that the poly3gauss basis allows you to add one extra basis function"
      write(*,*) "before similar instabilities develop."
      write(*,*) "-------------------------------------------------------------------------------"
      stop
    else if( alpha_max > 7 .and. print_message )then
      print_message = .false.
      write(*,*) "-------------------------------------------------------------------------------"
      write(*,*) "WARNING: Due to numerical instabilities in the basis construction for the poly3 <---- WARNING"
      write(*,*) "basis, it is strongly recommended not to exceed alpha_max = 7. For"
      write(*,*) "alpha_max > 10 this warning will turn into an error. You can do your own"
      write(*,*) "testing for 7 < alpha_max < 11, the results might still be useful within that"
      write(*,*) "range. Note that the poly3gauss basis allows you to add one extra basis function"
      write(*,*) "before similar instabilities develop."
      write(*,*) "-------------------------------------------------------------------------------"
    end if

    call radial_expansion_kernel(n_sites, n_neigh, rjs_in, alpha_max, rcut_soft_in, rcut_hard_in, &
                                 atom_sigma_in, atom_sigma_scaling, amplitude_scaling, nf, W, &
                                 scaling_mode, mask, radial_enhancement, do_derivatives, &
                                 do_central, central_weight, .false., exp_coeff, exp_coeff_der, &
                                 exp_coeff_der2)

  return
  end subroutine
!**************************************************************************




!**************************************************************************
!
! This subroutine returns the radial expansion coefficients using the
! polynomial basis set augmented with a Gaussian function at the origin
! and the saparable Gaussian representation for the atomic sites.
!
  subroutine get_radial_expansion_coefficients_poly3gauss(n_sites, n_neigh, rjs_in, alpha_max, rcut_soft_in, &
                                                          rcut_hard_in, atom_sigma_in, atom_sigma_scaling, &
                                                          amplitude_scaling, nf, W, scaling_mode, mask, &
                                                          radial_enhancement, do_derivatives, exp_coeff, &
                                                          exp_coeff_der, exp_coeff_der2)
!   Expansion coefficients using the polynomial basis with smooth filter plus a Gaussian centered at the origin

    implicit none

    integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
    real(dp), intent(in) :: rcut_soft_in, rcut_hard_in, rjs_in(:), atom_sigma_in, nf, atom_sigma_scaling
    real(dp), intent(in) :: amplitude_scaling
    logical, intent(in) :: mask(:), do_derivatives
    character(*), intent(in) :: scaling_mode
    real(dp) :: W(:,:)
    real(dp) :: exp_coeff(:,:), exp_coeff_der(:,:)
!   Optional second derivative d^2 exp_coeff / d rj^2, same shape. Requesting it
!   requires do_derivatives and legacy_filter_seed = .false.; see the kernel.
    real(dp), optional :: exp_coeff_der2(:,:)

    logical, save :: print_message = .true.

!   For this basis numerical instabilities in the orthonormal basis construction develop
!   above alpha_max = 8. These are relatively small up to alpha_max = 11 and become
!   catastrophic at alpha_max = 13.
    if( alpha_max > 11 )then
      write(*,*) "-------------------------------------------------------------------------------"
      write(*,*) "ERROR: Due to numerical instabilities in the basis construction for the         <---- ERROR"
      write(*,*) "poly3gauss basis, it is strongly recommended not to exceed alpha_max = 8. For"
      write(*,*) "alpha_max > 11 the instabilities are too large to proceed. You can do your own"
      write(*,*) "testing for 8 < alpha_max < 12, the results might still be useful within that"
      write(*,*) "range. Note that the poly3 basis allows you to add one *less* basis function"
      write(*,*) "before similar instabilities develop."
      write(*,*) "-------------------------------------------------------------------------------"
      stop
    else if( alpha_max > 8 .and. print_message )then
      print_message = .false.
      write(*,*) "-------------------------------------------------------------------------------"
      write(*,*) "WARNING: Due to numerical instabilities in the basis construction for the       <---- WARNING"
      write(*,*) "poly3gauss basis, it is strongly recommended not to exceed alpha_max = 8. For"
      write(*,*) "alpha_max > 11 this warning will turn into an error. You can do your own"
      write(*,*) "testing for 8 < alpha_max < 12, the results might still be useful within that"
      write(*,*) "range. Note that the poly3 basis allows you to add one *less* basis function"
      write(*,*) "before similar instabilities develop."
      write(*,*) "-------------------------------------------------------------------------------"
    end if

!   IMPORTANT: for this basis we skip i itself, which is neighbor number 1. That is what
!   do_central = .false. does in the kernel, so central_weight is never consulted.
    call radial_expansion_kernel(n_sites, n_neigh, rjs_in, alpha_max, rcut_soft_in, rcut_hard_in, &
                                 atom_sigma_in, atom_sigma_scaling, amplitude_scaling, nf, W, &
                                 scaling_mode, mask, radial_enhancement, do_derivatives, &
                                 .false., 0.d0, .true., exp_coeff, exp_coeff_der, &
                                 exp_coeff_der2)

  return
  end subroutine
!**************************************************************************




!**************************************************************************
!
! Shared kernel for the poly3 and poly3gauss radial expansions.
!
! Structure
! ---------
! The mathematics is unchanged from the original scalar version; what changes
! is the loop order.  The original walked the neighbour list one pair at a
! time and ran the alpha recursion inside, so the innermost loop carried a
! dependence and could not be vectorised, and every step called N_a() -- a
! square root -- afresh.
!
! Here the pair list is first compacted (cutoff and mask applied once) into
! blocks of BLK entries.  Every loop then runs over the block index, which is
! contiguous and dependence-free, with alpha as the *outer* loop.  The alpha
! recursion is still sequential in alpha, but each of its steps is now a
! vector operation over BLK pairs.  The N_a ratios that appear in the
! recursion and in the derivative are functions of alpha alone, so they are
! tabulated once per call rather than recomputed per pair.
!
! The two bases differ only in how many polynomial functions there are:
!   poly3       n_poly = alpha_max
!   poly3gauss  n_poly = alpha_max - 1, and slot alpha_max holds the overlap
!               with the Gaussian centred at the origin.
!
! Accuracy
! --------
! The filter recursion is seeded at alpha = -2, where (1 - r/rcut_hard)**(alpha+2)
! equals 1 rather than 0 at r = rcut_hard, so that first step carries a
! surface term the original dropped.  Keeping it makes exp_coeff_der agree
! with a finite difference of exp_coeff to machine precision; dropping it
! leaves a relative disagreement of order exp(-nf**2/2).  legacy_filter_seed
! selects which of the two is computed -- see the module header.
!
  subroutine radial_expansion_kernel(n_sites, n_neigh, rjs_in, alpha_max, rcut_soft_in, &
                                     rcut_hard_in, atom_sigma_in, atom_sigma_scaling, &
                                     amplitude_scaling, nf, W, scaling_mode, mask, &
                                     radial_enhancement, do_derivatives, do_central, &
                                     central_weight, with_gauss, exp_coeff, exp_coeff_der, &
                                     exp_coeff_der2)
    implicit none

    integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
    real(dp), intent(in) :: rcut_soft_in, rcut_hard_in, rjs_in(:), atom_sigma_in, nf, atom_sigma_scaling
    real(dp), intent(in) :: amplitude_scaling, central_weight
    logical, intent(in) :: mask(:), do_derivatives, do_central, with_gauss
    character(*), intent(in) :: scaling_mode
    real(dp), intent(in) :: W(:,:)
    real(dp), intent(inout) :: exp_coeff(:,:), exp_coeff_der(:,:)
    real(dp), intent(inout), optional :: exp_coeff_der2(:,:)

!   Pairs handled per vector block. The working set is
!   2*(alpha_max+6)*BLK + 5*alpha_max*BLK doubles, which for the alpha_max
!   values these bases support stays inside L1 at BLK = 64.
    integer, parameter :: BLK = 64

!   Buffers are automatic rather than allocatable so that BLK is the
!   compile-time leading dimension: T1(b,n) then costs one shift instead of a
!   multiply by a descriptor stride, and the b loops address contiguously.
    real(dp) :: T1(BLK, -3:alpha_max+4), T2(BLK, -3:alpha_max+4)
    real(dp) :: pol(BLK, alpha_max+4)
    real(dp) :: DT(BLK, alpha_max), dtot(BLK, alpha_max)
    real(dp) :: D2T(BLK, alpha_max), d2tot(BLK, alpha_max)
    real(dp) :: outv(BLK, alpha_max), outd(BLK, alpha_max), outh(BLK, alpha_max)
    real(dp) :: Nv(-3:alpha_max+4), ra(-1:alpha_max+4), rb(-1:alpha_max+4), rc(-1:alpha_max+4)
    real(dp) :: R1(alpha_max), R2(alpha_max), R3(alpha_max), R4(alpha_max)
    real(dp) :: Wl(alpha_max, alpha_max)

    real(dp) :: rj(BLK), sig(BLK), s2(BLK), isig(BLK), is2(BLK), is3(BLK)
    real(dp) :: amp(BLK), ampd(BLK), ampd2(BLK), acc(BLK), es(BLK), e0(BLK)
    real(dp) :: C1(BLK), C2(BLK), den(BLK), rcp(BLK)
    real(dp) :: pref(BLK), dpref(BLK), rjf(BLK), sgf(BLK), isgf(BLK), sf2(BLK), isf2(BLK)
    real(dp) :: iD(BLK), sqiD(BLK), drjf(BLK), dsgf(BLK)
    real(dp) :: ca(BLK), cb(BLK), cc(BLK)
    real(dp) :: g0(BLK), g1(BLK), g2(BLK)
    real(dp) :: B0(BLK), B1(BLK), B2(BLK), B3(BLK), B4(BLK)
    logical  :: is_central(BLK), keep(BLK)
    integer  :: idx(BLK)

    integer  :: n_poly, n_hi, i, j, k, n, nb
    real(dp) :: rcut_soft, atom_sigma, dr, inv_rcut, sqrt_rcut, inv_sqrt_rcut, inv_Nm2
    real(dp) :: pi, sq2, sqpi2, sq2pi, sq8pi, N_gauss, sigma_sq0, isq2
    real(dp) :: q, drnf
    logical  :: has_filter, do_poly, do_der2, unified

!   ------------------------------------------------------------------ setup
    pi    = dacos(-1.d0)
    sq2   = dsqrt(2.d0)
    isq2  = 1.d0/sq2
    sqpi2 = dsqrt(pi/2.d0)
    sq2pi = dsqrt(2.d0/pi)
    sq8pi = dsqrt(8.d0/pi)

    do_der2 = present(exp_coeff_der2)
    if( do_der2 .and. .not. do_derivatives )then
      write(*,*) "ERROR: soap_turbo_radial: exp_coeff_der2 was requested without do_derivatives"
      stop
    end if
!   The second derivative is built from the recurrence that the *exact*
!   integrals satisfy. With the legacy filter seed the stored integrals are not
!   those, so the result would be neither the derivative of what is returned
!   nor the derivative of anything else. Refuse rather than return a number
!   that looks plausible.
    if( do_der2 .and. legacy_filter_seed )then
      write(*,*) "-------------------------------------------------------------------------------"
      write(*,*) "ERROR: second radial derivatives need soap_turbo_radial's                <---- ERROR"
      write(*,*) "legacy_filter_seed = .false. . With the legacy seed the filter integrals are"
      write(*,*) "not the integrals the derivative recurrence assumes, so the second derivative"
      write(*,*) "would not be the second derivative of the coefficients that are returned."
      write(*,*) "-------------------------------------------------------------------------------"
      stop
    end if
!   With the corrected seed the whole coefficient -- plain Gaussian plus
!   filtered Gaussian -- obeys one derivative recurrence, because the filter is
!   a function of r alone and the rj dependence sits entirely in the Gaussian.
!   That is both simpler and cheaper than differentiating the two pieces
!   separately, so use it whenever it is available.
    unified = .not. legacy_filter_seed

    exp_coeff = 0.d0
    if( do_derivatives ) exp_coeff_der = 0.d0
    if( do_der2 ) exp_coeff_der2 = 0.d0

!   **************** New basis ******************
!   Redefine all the distances by dividing them by rcut_hard. We do this to
!   avoid numerical instability when the value of alpha is too high. From here
!   on rcut_hard == 1, which is why it no longer appears in the expressions.
    inv_rcut      = 1.d0 / rcut_hard_in
    rcut_soft     = rcut_soft_in * inv_rcut
    atom_sigma    = atom_sigma_in * inv_rcut
    dr            = 1.d0 - rcut_soft
    sqrt_rcut     = dsqrt(rcut_hard_in)
    inv_sqrt_rcut = 1.d0 / sqrt_rcut
!   *********************************************

    sigma_sq0 = atom_sigma*atom_sigma
    N_gauss   = dsqrt(2.d0/atom_sigma) / pi**0.25d0

!   A zero-width buffer (rcut_soft == rcut_hard, which read_files.f90 produces
!   from a buffer of 0) leaves the filter Gaussian with zero width. There is
!   nothing to integrate over [rcut_soft, rcut_hard] then, so skip the filter
!   rather than divide by dr -- the original divided, and returned NaN for
!   every coefficient.
    has_filter = ( dr > 0.d0 )
    q    = 0.d0
    drnf = 0.d0
    if( has_filter )then
      drnf = dr / nf
      q    = drnf * drnf
    end if

    if( with_gauss )then
      n_poly = alpha_max - 1
    else
      n_poly = alpha_max
    end if
    do_poly = ( n_poly >= 1 )
    n_hi = n_poly
    if( do_derivatives ) n_hi = n_poly + 2
    if( do_der2 )        n_hi = n_poly + 4

!   W arrives as a slice of a larger, species-blocked matrix, so it is not
!   contiguous. Copy it once; it is only alpha_max**2 elements.
    Wl(1:alpha_max, 1:alpha_max) = W(1:alpha_max, 1:alpha_max)

!   ------------------------------------------------- alpha-only coefficients
!   N_a(rcut_hard=1, a) = sqrt(1/(2a+5)). The original called N_a() -- a
!   square root -- once per alpha per pair in the recursion and four times per
!   alpha per pair in the derivative; these are functions of alpha alone.
!   Index -3 is the dead slot that the recursion multiplies by (n+1) = 0 on
!   its first step; the original seeded it with 1.
    Nv(-3) = 1.d0
    do n = -2, alpha_max+4
      Nv(n) = dsqrt( 1.d0 / dfloat(2*n+5) )
    end do
    do n = -1, n_hi
      ra(n) = dfloat(n+1) * Nv(n-2) / Nv(n)
      rb(n) = Nv(n-1) / Nv(n)
      rc(n) = 1.d0 / Nv(n)
    end do
    ra(-1)  = 0.d0
    inv_Nm2 = 1.d0 / Nv(-2)

    if( do_derivatives .and. do_poly )then
      do n = 1, n_poly
        R1(n) = Nv(n+1) / Nv(n)
        R2(n) = Nv(n+2) / Nv(n)
      end do
      if( do_der2 )then
        do n = 1, n_poly
          R3(n) = Nv(n+3) / Nv(n)
          R4(n) = Nv(n+4) / Nv(n)
        end do
      end if
    end if

!   ------------------------------------------------- walk the neighbour list
    nb = 0
    k  = 0
    do i = 1, n_sites
      do j = 1, n_neigh(i)
        k = k + 1
        if( j == 1 .and. .not. do_central ) cycle
        if( rjs_in(k) >= rcut_hard_in ) cycle
        if( .not. mask(k) ) cycle
        nb = nb + 1
        idx(nb)        = k
        rj(nb)         = rjs_in(k) * inv_rcut
        is_central(nb) = ( j == 1 )
        if( nb == BLK )then
          call process_block()
          nb = 0
        end if
      end do
    end do
    if( nb > 0 ) call process_block()

  return

  contains

!   Everything below runs on the nb pairs currently gathered in idx(1:nb).
!   Every loop over b is contiguous and free of loop-carried dependences; the
!   alpha recursion is sequential in alpha but vector over the block.
    subroutine process_block()
      integer :: b, n, m, kk
      real(dp) :: f, g, h, hp, hpp, wnm, ss2, iss, siss, ex2, aa
      real(dp) :: is4, is5, is6, k0, k1, k2, k3, k4
      real(dp) :: M0, M1, M2, M3, M4, pp, sc0, ge
      logical  :: gm

!     ------------------------- per-pair scalars; reciprocals in their own
!     loop so that the divides vectorise instead of sitting in a dependence
!     chain with the multiplies that use them
      do b = 1, nb
        sig(b) = atom_sigma + atom_sigma_scaling*rj(b)
        s2(b)  = sig(b)*sig(b)
      end do
      do b = 1, nb
        isig(b) = 1.d0/sig(b)
      end do
      do b = 1, nb
        is2(b) = isig(b)*isig(b)
        is3(b) = is2(b)*isig(b)
      end do

!     ------------------------------------------------------------ amplitude
!     WARNING: the 1/atom_sigma_angular^2 term is missing from these amplitudes
!     and needs to be taken into account in the corresponding part of the code.
!
!     WARNING2: these expressions already assume rcut_hard = 1, so that
!     parameter is missing from them.
!
!     amplitude_scaling and scaling_mode are per-call scalars, so the branch on
!     them is hoisted out of the block loop; the per-pair test on the cubic is a
!     mask rather than a branch, to keep the loop vectorisable.
!
!     ampd2 is only filled when second derivatives are wanted. With
!     h(rj) = 1 + 2rj^3 - 3rj^2 and sigma = sigma_0 + c*rj,
!       a    = h^s / sigma
!       a'   = s h^(s-1) h' / sigma - c h^s / sigma^2
!       a''  = [s(s-1) h^(s-2) h'^2 + s h^(s-1) h''] / sigma
!              - 2c s h^(s-1) h' / sigma^2 + 2 c^2 h^s / sigma^3
      if( scaling_mode /= "polynomial" )then
!       The original left amplitude and amplitude_der undefined here.
!       read_files.f90 forces scaling_mode = "polynomial", so this branch is
!       unreachable from an input file; an unscaled amplitude is at least a
!       defined answer.
        do b = 1, nb
          amp(b)   = 1.d0
          ampd(b)  = 0.d0
          ampd2(b) = 0.d0
        end do
      else if( amplitude_scaling == 0.d0 )then
        do b = 1, nb
          amp(b)   = isig(b)
          ampd(b)  = -atom_sigma_scaling*is2(b)
          ampd2(b) = 2.d0*atom_sigma_scaling*atom_sigma_scaling*is3(b)
        end do
      else if( amplitude_scaling == 1.d0 )then
        do b = 1, nb
          h        = 1.d0 + 2.d0*rj(b)**3 - 3.d0*rj(b)**2
          hp       = 6.d0*rj(b)**2 - 6.d0*rj(b)
          hpp      = 12.d0*rj(b) - 6.d0
          g        = isig(b)*h
          amp(b)   = merge( 0.d0, g, h <= 1.d-10 )
          ampd(b)  = merge( 0.d0, hp*isig(b) - atom_sigma_scaling*isig(b)*g, h <= 1.d-10 )
          ampd2(b) = merge( 0.d0, hpp*isig(b) - 2.d0*atom_sigma_scaling*hp*is2(b) &
                                  + 2.d0*atom_sigma_scaling*atom_sigma_scaling*h*is3(b), &
                            h <= 1.d-10 )
        end do
      else
        do b = 1, nb
          h = 1.d0 + 2.d0*rj(b)**3 - 3.d0*rj(b)**2
          if( h <= 1.d-10 )then
            amp(b)   = 0.d0
            ampd(b)  = 0.d0
            ampd2(b) = 0.d0
          else
            hp       = 6.d0*rj(b)**2 - 6.d0*rj(b)
            hpp      = 12.d0*rj(b) - 6.d0
            aa       = amplitude_scaling
            amp(b)   = isig(b) * h**aa
            ampd(b)  = aa*h**(aa-1.d0)*hp*isig(b) - atom_sigma_scaling*isig(b)*amp(b)
            ampd2(b) = ( aa*(aa-1.d0)*h**(aa-2.d0)*hp*hp + aa*h**(aa-1.d0)*hpp )*isig(b) &
                       - 2.d0*atom_sigma_scaling*aa*h**(aa-1.d0)*hp*is2(b) &
                       + 2.d0*atom_sigma_scaling*atom_sigma_scaling*h**aa*is3(b)
          end if
        end do
      end if

!     The central atom needs to be scaled by central_weight
      if( do_central .and. central_weight /= 1.d0 )then
        do b = 1, nb
          if( is_central(b) )then
            amp(b)   = central_weight * amp(b)
            ampd(b)  = central_weight * ampd(b)
            ampd2(b) = central_weight * ampd2(b)
          end if
        end do
      end if

!     The radial enhancement adds a scaling corresponding to the integral of a
!     Gaussian at the position of atom j. It multiplies the amplitude by a
!     factor m(rj), so (a m)'' = a'' m + 2 a' m' + a m''.
      if( radial_enhancement == 1 )then
        do b = 1, nb
          g        = rj(b) + sq2pi*sig(b)
          f        = 1.d0 + sq2pi*atom_sigma_scaling
          ampd2(b) = ampd2(b)*g + 2.d0*ampd(b)*f
          ampd(b)  = amp(b)*f + ampd(b)*g
          amp(b)   = amp(b)*g
        end do
      else if( radial_enhancement == 2 )then
        do b = 1, nb
          g        = rj(b)**2 + s2(b) + sq8pi*sig(b)*rj(b)
          f        = 2.d0*rj(b) + 2.d0*sig(b)*atom_sigma_scaling &
                     + sq8pi*sig(b) + sq8pi*rj(b)*atom_sigma_scaling
          h        = 2.d0 + 2.d0*atom_sigma_scaling*atom_sigma_scaling &
                     + 2.d0*sq8pi*atom_sigma_scaling
          ampd2(b) = ampd2(b)*g + 2.d0*ampd(b)*f + amp(b)*h
          ampd(b)  = amp(b)*f + ampd(b)*g
          amp(b)   = amp(b)*g
        end do
      end if

!     ------------------------- overlap with the Gaussian on [0, rcut_soft]
!     C1 is the surface term at rcut_soft, C2 the one at the origin; both pick
!     up one more power of their cutoff per step. es and e0 are kept because the
!     Gaussian-at-the-origin block below needs e0 again.
      do b = 1, nb
        es(b) = dexp( -0.5d0*(rcut_soft-rj(b))**2 * is2(b) )
        e0(b) = dexp( -0.5d0*rj(b)**2 * is2(b) )
      end do
      do b = 1, nb
        C1(b)    = s2(b)*es(b)
        C2(b)    = s2(b)*e0(b)
        T1(b,-3) = 0.d0
        T1(b,-2) = sqpi2 * sig(b) * ( derf( (rcut_soft-rj(b))*isq2*isig(b) ) &
                                    + derf( rj(b)*isq2*isig(b) ) ) * inv_Nm2
      end do

      do n = -1, n_hi
        if( n >= 0 )then
          do b = 1, nb
            C1(b) = C1(b)*dr
          end do
        end if
        do b = 1, nb
          T1(b,n) = s2(b)*ra(n)*T1(b,n-2) - rb(n)*(rj(b)-1.d0)*T1(b,n-1) &
                    + ( C1(b) - C2(b) )*rc(n)
        end do
      end do

!     ------- legacy derivative of that piece alone (amplitude excluded)
      if( do_derivatives .and. do_poly .and. .not. unified )then
        do b = 1, nb
          f     = rj(b) - 1.d0
          g     = atom_sigma_scaling * f * isig(b)
          ca(b) = f*is2(b)*( g - 1.d0 )
          cb(b) = is2(b)*( 2.d0*g - 1.d0 )
          cc(b) = atom_sigma_scaling*is3(b)
        end do
        do n = 1, n_poly
          do b = 1, nb
            DT(b,n) = ca(b)*T1(b,n) + R1(n)*cb(b)*T1(b,n+1) + R2(n)*cc(b)*T1(b,n+2)
          end do
        end do
      end if

!     ------------------------- Gaussian filter on [rcut_soft, rcut_hard]
!     If the atom is less than 4 standard deviations away from the soft cutoff
!     we add this correction to the integrals. The product of the atomic
!     Gaussian and the filter is another Gaussian, of width sgf about rjf, times
!     the prefactor pref.
      if( has_filter )then
        do b = 1, nb
          den(b) = s2(b) + q
        end do
        do b = 1, nb
          iD(b) = 1.d0/den(b)
        end do
        do b = 1, nb
          sqiD(b) = dsqrt(iD(b))
          sgf(b)  = sig(b)*drnf*sqiD(b)
          rjf(b)  = ( s2(b)*rcut_soft + q*rj(b) )*iD(b)
          sf2(b)  = sgf(b)*sgf(b)
!         outside 4 sigma the original drops the whole filter contribution
          keep(b) = ( (rcut_soft - rj(b)) < 4.d0*sig(b) )
        end do
        do b = 1, nb
          isgf(b) = 1.d0/sgf(b)
        end do
        do b = 1, nb
          isf2(b) = isgf(b)*isgf(b)
          pref(b) = merge( dexp( -0.5d0*(rcut_soft-rj(b))**2 * iD(b) ), 0.d0, keep(b) )
        end do

        do b = 1, nb
          C2(b)    = sf2(b) * dexp( -0.5d0*(rcut_soft-rjf(b))**2 * isf2(b) )
          T2(b,-3) = 0.d0
          T2(b,-2) = sqpi2 * sgf(b) * ( derf( (1.d0-rjf(b))*isq2*isgf(b) ) &
                                      - derf( (rcut_soft-rjf(b))*isq2*isgf(b) ) ) * inv_Nm2
        end do

        do n = -1, n_hi
          if( n >= 0 )then
            do b = 1, nb
              C2(b) = C2(b)*dr
            end do
          end if
          do b = 1, nb
            T2(b,n) = sf2(b)*ra(n)*T2(b,n-2) - rb(n)*(rjf(b)-1.d0)*T2(b,n-1) - C2(b)*rc(n)
          end do
          if( n == -1 .and. unified )then
!           Surface term at r = rcut_hard. (1-r/rcut_hard)**(alpha+2) is zero
!           there for every alpha >= -1, but this step is alpha = -2, where the
!           exponent is 0 and the factor is 1. Dropping it is what makes the
!           original's exp_coeff_der disagree with a central difference of its
!           own exp_coeff.
            do b = 1, nb
              T2(b,-1) = T2(b,-1) + sf2(b)*dexp( -0.5d0*(1.d0-rjf(b))**2 * isf2(b) )*rc(-1)
            end do
          end if
        end do

!       ------- legacy derivative of the filter piece, differentiated through
!       pref, rj_f and sigma_f separately
        if( do_derivatives .and. do_poly .and. .not. unified )then
          do b = 1, nb
            f        = rcut_soft - rj(b)
            dpref(b) = pref(b) * ( f*iD(b) + f*f*iD(b)*iD(b)*sig(b)*atom_sigma_scaling )
            drjf(b)  = ( 2.d0*sig(b)*rcut_soft*atom_sigma_scaling + q )*iD(b) &
                       - ( s2(b)*rcut_soft + q*rj(b) )*2.d0*sig(b)*atom_sigma_scaling*iD(b)*iD(b)
            dsgf(b)  = atom_sigma_scaling*drnf*sqiD(b)*( 1.d0 - s2(b)*iD(b) )
          end do
          do b = 1, nb
            f     = rjf(b) - 1.d0
            g     = dsgf(b) * f * isgf(b)
            ca(b) = pref(b) * f*isf2(b)*( g - drjf(b) )
            cb(b) = pref(b) * isf2(b)*( 2.d0*g - drjf(b) )
            cc(b) = pref(b) * dsgf(b)*isgf(b)*isf2(b)
          end do
          do n = 1, n_poly
            do b = 1, nb
              DT(b,n) = DT(b,n) + ca(b)*T2(b,n) + R1(n)*cb(b)*T2(b,n+1) &
                        + R2(n)*cc(b)*T2(b,n+2) + dpref(b)*T2(b,n)
            end do
          end do
        end if
      else
        do b = 1, nb
          pref(b) = 0.d0
        end do
        do n = 1, n_hi
          do b = 1, nb
            T2(b,n) = 0.d0
          end do
        end do
      end if

!     --------------- the polynomial part of the coefficient, both pieces
      do n = 1, n_hi
        do b = 1, nb
          pol(b,n) = T1(b,n) + pref(b)*T2(b,n)
        end do
      end do

!     --------------- unified derivatives
!     The filter is a function of r alone, so the whole rj dependence of
!     C_alpha = int g_alpha(r) exp(-(r-rj)^2/2 sigma(rj)^2) F(r) dr sits in the
!     Gaussian.  With u = r - rj,
!         dG/drj   = [ u/sigma^2 + c u^2/sigma^3 ] G
!         d2G/drj2 = [ -1/sigma^2 - 4c u/sigma^3 + (1-3c^2) u^2/sigma^4
!                      + 2c u^3/sigma^5 + c^2 u^4/sigma^6 ] G
!     and u^m is turned into a combination of C_(alpha+p) by
!         (rcut_hard - r)^p g_alpha = rcut_hard^p (N_(alpha+p)/N_alpha) g_(alpha+p),
!     writing u = (rcut_hard - rj) - (rcut_hard - r).  So one recurrence covers
!     both pieces, and the second derivative reaches alpha+4.
      if( do_derivatives .and. do_poly .and. unified )then
        do b = 1, nb
          f     = 1.d0 - rj(b)
          g0(b) = f*is2(b) + atom_sigma_scaling*f*f*is3(b)
          g1(b) = -( is2(b) + 2.d0*atom_sigma_scaling*f*is3(b) )
          g2(b) = atom_sigma_scaling*is3(b)
        end do
        do n = 1, n_poly
          do b = 1, nb
            DT(b,n) = g0(b)*pol(b,n) + R1(n)*g1(b)*pol(b,n+1) + R2(n)*g2(b)*pol(b,n+2)
          end do
        end do

        if( do_der2 )then
          do b = 1, nb
            f   = 1.d0 - rj(b)
            is4 = is2(b)*is2(b)
            is5 = is4*isig(b)
            is6 = is4*is2(b)
            k0  = -is2(b)
            k1  = -4.d0*atom_sigma_scaling*is3(b)
            k2  = ( 1.d0 - 3.d0*atom_sigma_scaling*atom_sigma_scaling )*is4
            k3  = 2.d0*atom_sigma_scaling*is5
            k4  = atom_sigma_scaling*atom_sigma_scaling*is6
            B0(b) = k0 + k1*f + k2*f*f + k3*f**3 + k4*f**4
            B1(b) = -( k1 + 2.d0*k2*f + 3.d0*k3*f*f + 4.d0*k4*f**3 )
            B2(b) = k2 + 3.d0*k3*f + 6.d0*k4*f*f
            B3(b) = -( k3 + 4.d0*k4*f )
            B4(b) = k4
          end do
          do n = 1, n_poly
            do b = 1, nb
              D2T(b,n) = B0(b)*pol(b,n) + R1(n)*B1(b)*pol(b,n+1) + R2(n)*B2(b)*pol(b,n+2) &
                         + R3(n)*B3(b)*pol(b,n+3) + R4(n)*B4(b)*pol(b,n+4)
            end do
          end do
        end if
      end if

!     ------------------------- overlap with the Gaussian at the origin
!     Only for poly3gauss, and only if the atom is close enough to the origin
!     that the overlap is not negligible. The filter is not applied to it, which
!     is a good approximation as long as rcut_soft >= 4*atom_sigma.
!
!     This one is not a polynomial, so the alpha-shift trick above does not
!     apply. Instead the moments M_m = int phi(r) (r-rj)^m G dr obey
!       M_m = [ sigma_0^2 sigma^2 N_gauss (-rj)^(m-1) e0
!               - sigma^2 rj M_(m-1) + sigma_0^2 sigma^2 (m-1) M_(m-2) ] / sigma_star^2
!     seeded by M_0, the overlap itself, and the derivatives follow from the
!     same dG/drj and d2G/drj2 as above.
      if( with_gauss )then
        do b = 1, nb
          ss2  = sigma_sq0 + s2(b)
          iss  = 1.d0/ss2
          siss = dsqrt(iss)
!         beyond 4 sigma the original drops the overlap; zeroing the surface
!         term as well keeps every moment, hence every derivative, zero too
          gm   = ( rj(b) < 4.d0*(atom_sigma + sig(b)) )
          ge   = merge( N_gauss*e0(b), 0.d0, gm )
          M0   = merge( dexp( -0.5d0*rj(b)**2*iss ) * sqpi2 * sig(b)*atom_sigma*siss &
                        * ( 1.d0 + derf( atom_sigma*isig(b)*rj(b)*isq2*siss ) ) * N_gauss, &
                        0.d0, gm )
          pol(b,alpha_max) = M0
          if( do_derivatives )then
            sc0 = sigma_sq0*s2(b)*iss
            pp  = s2(b)*rj(b)*iss
!           e0 is exp(-rj^2/(2 sigma^2)), which is what phi(0)*G(0) needs
            M1 = sc0*ge - pp*M0
            M2 = sc0*ge*(-rj(b)) - pp*M1 + sc0*M0
            DT(b,alpha_max) = M1*is2(b) + atom_sigma_scaling*M2*is3(b)
            if( do_der2 )then
              M3  = sc0*ge*rj(b)**2 - pp*M2 + 2.d0*sc0*M1
              M4  = sc0*ge*(-rj(b)**3) - pp*M3 + 3.d0*sc0*M2
              is4 = is2(b)*is2(b)
              is5 = is4*isig(b)
              is6 = is4*is2(b)
              D2T(b,alpha_max) = -M0*is2(b) - 4.d0*atom_sigma_scaling*M1*is3(b) &
                   + ( 1.d0 - 3.d0*atom_sigma_scaling*atom_sigma_scaling )*M2*is4 &
                   + 2.d0*atom_sigma_scaling*M3*is5 &
                   + atom_sigma_scaling*atom_sigma_scaling*M4*is6
            end if
          end if
        end do
      end if

!     ------------------------ transform from g_alpha to g_n and store
      if( do_derivatives )then
        do n = 1, alpha_max
          do b = 1, nb
            dtot(b,n) = amp(b)*DT(b,n) + ampd(b)*pol(b,n)
          end do
        end do
        if( do_der2 )then
          do n = 1, alpha_max
            do b = 1, nb
              d2tot(b,n) = amp(b)*D2T(b,n) + 2.d0*ampd(b)*DT(b,n) + ampd2(b)*pol(b,n)
            end do
          end do
        end if
        do n = 1, alpha_max
          wnm = Wl(n,1)
          do b = 1, nb
            acc(b) = wnm*dtot(b,1)
          end do
          do m = 2, alpha_max
            wnm = Wl(n,m)
            do b = 1, nb
              acc(b) = acc(b) + wnm*dtot(b,m)
            end do
          end do
!         The original applied the 1/sqrt(rcut_hard) from the change of variable
!         in a separate sweep over the whole array afterwards; folding it in
!         here saves a read-modify-write pass over alpha_max*n_atom_pairs
!         doubles.
          do b = 1, nb
            outd(b,n) = acc(b)*inv_sqrt_rcut
          end do
        end do
        if( do_der2 )then
          do n = 1, alpha_max
            wnm = Wl(n,1)
            do b = 1, nb
              acc(b) = wnm*d2tot(b,1)
            end do
            do m = 2, alpha_max
              wnm = Wl(n,m)
              do b = 1, nb
                acc(b) = acc(b) + wnm*d2tot(b,m)
              end do
            end do
            do b = 1, nb
              outh(b,n) = acc(b)*inv_sqrt_rcut*inv_rcut
            end do
          end do
        end if
      end if

      do n = 1, alpha_max
        wnm = Wl(n,1)
        do b = 1, nb
          acc(b) = wnm*pol(b,1)
        end do
        do m = 2, alpha_max
          wnm = Wl(n,m)
          do b = 1, nb
            acc(b) = acc(b) + wnm*pol(b,m)
          end do
        end do
        do b = 1, nb
          outv(b,n) = amp(b)*acc(b)*sqrt_rcut
        end do
      end do

!     exp_coeff(:,k) is one contiguous column per pair
      do b = 1, nb
        kk = idx(b)
        do n = 1, alpha_max
          exp_coeff(n,kk) = outv(b,n)
        end do
      end do
      if( do_derivatives )then
        do b = 1, nb
          kk = idx(b)
          do n = 1, alpha_max
            exp_coeff_der(n,kk) = outd(b,n)
          end do
        end do
      end if
      if( do_der2 )then
        do b = 1, nb
          kk = idx(b)
          do n = 1, alpha_max
            exp_coeff_der2(n,kk) = outh(b,n)
          end do
        end do
      end if

    end subroutine process_block

  end subroutine radial_expansion_kernel
!**************************************************************************










!**************************************************************************
!
! This subroutine returns the overlap matrix S and the orthonormalization
! matrix W = S^{-1/2} needed to construct the orthonormal basis from the
! polynomial basis set. It is also needed to transform between original
! basis and orthonormal basis. It requires blas/lapack to work.
!
  subroutine get_orthonormalization_matrix_poly3(alpha_max, S, W)

    implicit none

    integer :: alpha_max, i, j, info
    real(dp), intent(inout) :: W(:,:), S(:,:)
    real(dp), allocatable :: Sb(:,:), U(:,:), VT(:,:), svd(:), work(:), Sc(:,:)
    integer, allocatable :: ipiv(:)
    logical :: stable_basis = .true.

    allocate( Sb(1:alpha_max, 1:alpha_max) )
    allocate( Sc(1:alpha_max, 1:alpha_max) )
    allocate( U(1:alpha_max, 1:alpha_max) )
    allocate( Vt(1:alpha_max, 1:alpha_max) )
    allocate( svd(1:alpha_max) )
    allocate( work(1:6*alpha_max) )
    allocate( ipiv(1:alpha_max) )

    do i = 1, alpha_max
      S(i,i) = 1.d0
      do j = i+1, alpha_max
        S(i,j) = dsqrt( dfloat(5+2*i) * dfloat(5+2*j) ) / dfloat(5+i+j)
        S(j,i) = S(i,j)
      end do
    end do

    Sb(1:alpha_max, 1:alpha_max) = S(1:alpha_max, 1:alpha_max)

!   Do Singular Value Decomposition of S
    call dgesvd( "A", "A", alpha_max, alpha_max, S, alpha_max, svd, U, alpha_max, VT, &
                alpha_max, work, 6*alpha_max, info )
!   For debugging
    if( .false. )then
      do i = 1, alpha_max
        write(*,*) i, svd(i)
      end do
    end if
!   S^0.5
    S = 0.d0
    do i = 1, alpha_max
      S(i,i) = dsqrt(svd(i))
    end do
!    S = matmul(U,S)
!    S = matmul(S,VT)
    call dgemm("N", "N", alpha_max, alpha_max, alpha_max, 1.d0, U, alpha_max, S, alpha_max, 0.d0, Sc, alpha_max)
    call dgemm("N", "N", alpha_max, alpha_max, alpha_max, 1.d0, Sc, alpha_max, VT, alpha_max, 0.d0, S, alpha_max)
!   Invert S
    if( stable_basis )then
      call dpotrf( "U", alpha_max, S, alpha_max, info )
      call dpotri( "U", alpha_max, S, alpha_max, info )
    else
!     These ones are very unstable
!      call dgetrf( alpha_max, alpha_max, S, alpha_max, ipiv, info )
!      call dgetri( alpha_max, S, alpha_max, ipiv, work, 6*alpha_max, info )
!     These are essentially the same as dpotr*
      call dsytrf( "U", alpha_max, S, alpha_max, ipiv, work, 6*alpha_max, info )
      call dsytri( "U", alpha_max, S, alpha_max, ipiv, work(1:alpha_max), info )
    end if
    do i = 1, alpha_max
      W(i,i) = S(i,i)
      do j = i+1, alpha_max
        W(i,j) = S(i,j)
        W(j,i) = S(i,j)
      end do
    end do

    S(1:alpha_max, 1:alpha_max) = Sb(1:alpha_max, 1:alpha_max)

    deallocate( Sb, U, Vt, svd, work, ipiv, Sc )

  return
  end subroutine
!**************************************************************************









!**************************************************************************
!
! This subroutine returns the overlap matrix S and the orthonormalization
! matrix W = S^{-1/2} needed to construct the orthonormal basis from the
! polynomial basis set augmented with a central Gaussian function. It is
! also needed to transform between original basis and orthonormal basis.
! It requires blas/lapack to work.
!
  subroutine get_orthonormalization_matrix_poly3gauss(alpha_max, atom_sigma_in, rcut_hard_in, S, W)

    implicit none

    integer :: alpha_max, i, j, info, n
    real(dp), intent(inout) :: W(:,:), S(:,:)
    real(dp), allocatable :: Sb(:,:), U(:,:), VT(:,:), svd(:), work(:), Sc(:,:)
    real(dp) :: s2, I_n, N_n, N_np1, I_np1, N_np2, I_np2, C2, sq2, pi, atom_sigma, rcut_hard
    real(dp), intent(in) :: rcut_hard_in, atom_sigma_in
    integer, allocatable :: ipiv(:)
    logical :: stable_basis = .true.

    allocate( Sb(1:alpha_max, 1:alpha_max) )
    allocate( Sc(1:alpha_max, 1:alpha_max) )
    allocate( U(1:alpha_max, 1:alpha_max) )
    allocate( Vt(1:alpha_max, 1:alpha_max) )
    allocate( svd(1:alpha_max) )
    allocate( work(1:6*alpha_max) )
    allocate( ipiv(1:alpha_max) )

!   These are the overlap integrals for the polynomial functions
    do i = 1, alpha_max-1
      S(i,i) = 1.d0
      do j = i+1, alpha_max-1
        S(i,j) = dsqrt( dfloat(5+2*i) * dfloat(5+2*j) ) / dfloat(5+i+j)
        S(j,i) = S(i,j)
      end do
    end do

!   These are the overlap integrals between the Gaussian and the polynomials
!   See derivation of radial expansion coefficients to understand this code
!   **** New basis ****
    atom_sigma = atom_sigma_in/rcut_hard_in
    rcut_hard = 1.d0
!   *******************
    s2 = atom_sigma**2
    sq2 = dsqrt(2.d0)
    pi = dacos(-1.d0)
    I_n = 0.d0
    N_n = 1.d0
    N_np1 = N_a(rcut_hard, -2)
    I_np1 = dsqrt(pi/2.d0) * atom_sigma * derf( rcut_hard/sq2/atom_sigma ) / N_np1
    C2 = s2 / rcut_hard
    do n = -1, alpha_max-1
      C2 = C2 * rcut_hard
      N_np2 = N_a(rcut_hard, n)
      I_np2 = s2 * dfloat(n+1) * N_n/ N_np2 * I_n &
              + N_np1 * rcut_hard / N_np2 * I_np1 &
              - C2 / N_np2
      if(n > 0)then
!       Include the normalization factor of the Gaussian
        S(alpha_max, n) = I_np2 * sq2 / dsqrt(atom_sigma) / pi**0.25d0
        S(n, alpha_max) = S(alpha_max, n)
      end if
      N_n = N_np1
      N_np1 = N_np2
      I_n = I_np1
      I_np1 = I_np2
    end do
    S(alpha_max, alpha_max) = 1.d0

    Sb(1:alpha_max, 1:alpha_max) = S(1:alpha_max, 1:alpha_max)

!   Do Singular Value Decomposition of S
    call dgesvd( "A", "A", alpha_max, alpha_max, S, alpha_max, svd, U, alpha_max, VT, &
                alpha_max, work, 6*alpha_max, info )
!   For debugging
    if( .false. )then
      do i = 1, alpha_max
        write(*,*) i, svd(i)
      end do
    end if
!   S^0.5
    S = 0.d0
    do i = 1, alpha_max
      S(i,i) = dsqrt(svd(i))
    end do
!    S = matmul(U,S)
!    S = matmul(S,VT)
    call dgemm("N", "N", alpha_max, alpha_max, alpha_max, 1.d0, U, alpha_max, S, alpha_max, 0.d0, Sc, alpha_max)
    call dgemm("N", "N", alpha_max, alpha_max, alpha_max, 1.d0, Sc, alpha_max, VT, alpha_max, 0.d0, S, alpha_max)
!   Invert S
    if( stable_basis )then
      call dpotrf( "U", alpha_max, S, alpha_max, info )
      call dpotri( "U", alpha_max, S, alpha_max, info )
    else
!     These ones are very unstable
!      call dgetrf( alpha_max, alpha_max, S, alpha_max, ipiv, info )
!      call dgetri( alpha_max, S, alpha_max, ipiv, work, 6*alpha_max, info )
!     These are essentially the same as dpotr*
      call dsytrf( "U", alpha_max, S, alpha_max, ipiv, work, 6*alpha_max, info )
      call dsytri( "U", alpha_max, S, alpha_max, ipiv, work(1:alpha_max), info )
    end if
    do i = 1, alpha_max
      W(i,i) = S(i,i)
      do j = i+1, alpha_max
        W(i,j) = S(i,j)
        W(j,i) = S(i,j)
      end do
    end do

    S(1:alpha_max, 1:alpha_max) = Sb(1:alpha_max, 1:alpha_max)

    deallocate( Sb, U, Vt, svd, work, ipiv, Sc )

  return
  end subroutine
!**************************************************************************






!**************************************************************************
!
! This subroutine returns the overlap matrix S and the orthonormalization
! matrix W = S^{-1/2} needed to construct the orthonormal basis from the
! polynomial basis set. It is also needed to transform between original
! basis and orthonormal basis. Unlike the original function, it does not
! require blas/lapack to work and instead relies on pretabulated values.
!
  subroutine get_orthonormalization_matrix_poly3_tabulated(alpha_max, S, W)

    implicit none

    integer :: alpha_max
    real(dp), intent(inout) :: W(:,:), S(:,:)

    select case(alpha_max)
      case(:0)
        write(*,*) "Bad value of alpha_max"
        stop
      case(1)
        W = reshape([1.0000000000000000d0], shape(W))
        S = reshape([1.0000000000000000d0], shape(S))
      case(2)
        W = reshape([5.9999999999999920d0, -5.2915026221291726d0, -5.2915026221291726d0,  &
                    5.9999999999999920d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.99215674164922152d0,  &
                    1.0000000000000000d0], shape(S))
      case(3)
        W = reshape([16.532871534059733d0, -29.078310649176263d0, 13.308493852760057d0, -29.078310649176263d0,  &
                    65.256761373659685d0, -36.000096455630107d0, 13.308493852760057d0, -36.000096455630107d0,  &
                    23.492063480194044d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0], shape(S))
      case(4)
        W = reshape([31.619473929163473d0, -82.899603631697133d0, 76.876069386827055d0, -24.858289195076921d0,  &
                    -82.899603631697133d0, 278.76362329353776d0, -309.81360154834925d0, 114.16667789081102d0,  &
                    76.876069386827055d0, -309.81360154834925d0, 401.66898510515369d0, -168.42692378140524d0,  &
                    -24.858289195076921d0, 114.16667789081102d0, -168.42692378140524d0, 79.877446522677999d0 &
                    ], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0,  &
                    0.98333216603563356d0, 0.97499604304356913d0, 0.99498743710661997d0, 1.0000000000000000d0,  &
                    0.99652172859178323d0, 0.95393920141694566d0, 0.98333216603563356d0, 0.99652172859178323d0,  &
                    1.0000000000000000d0], shape(S))
      case(5)
        W = reshape([48.701840798683406d0, -167.01751270577688d0, 232.29114353203897d0, -152.17733678177117d0,  &
                    38.937949480633065d0, -167.01751270577688d0, 745.79170301284648d0, -1251.5498958156577d0,  &
                    936.75356599200450d0, -263.84749322361466d0, 232.29114353203897d0, -1251.5498958156577d0,  &
                    2447.7952379389576d0, -2072.0455670547381d0, 643.96375585768624d0, -152.17733678177117d0,  &
                    936.75356599200450d0, -2072.0455670547381d0, 1959.7497439691933d0, -672.11823364760107d0,  &
                    38.937949480633065d0, -263.84749322361466d0, 643.96375585768624d0, -672.11823364760107d0,  &
                    253.84463194408610d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.99215674164922152d0, 1.0000000000000000d0,  &
                    0.99498743710661997d0, 0.98333216603563356d0, 0.96824583655185414d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0, 0.99652172859178323d0, 0.98809481374347152d0,  &
                    0.95393920141694566d0, 0.98333216603563356d0, 0.99652172859178323d0, 1.0000000000000000d0,  &
                    0.99744571741206722d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0], shape(S))
      case(6)
        W = reshape([65.344980865540876d0, -270.96483658090318d0, 495.43665447034391d0, -487.28339480337212d0,  &
                    252.44181841585481d0, -54.246284684757562d0, -270.96483658090318d0, 1493.2157334161463d0,  &
                    -3338.9403306650820d0, 3783.7603589943619d0, -2167.7466416209732d0, 500.79512237094178d0,  &
                    495.43665447034391d0, -3338.9403306650820d0, 8760.0319024860873d0, -11249.630736045141d0,  &
                    7104.9320746203475d0, -1771.4461242822615d0, -487.28339480337212d0, 3783.7603589943619d0,  &
                    -11249.630736045141d0, 16097.101050349149d0, -11150.857982565711d0, 3007.1993994548748d0,  &
                    252.44181841585481d0, -2167.7466416209732d0, 7104.9320746203475d0, -11150.857982565711d0,  &
                    8419.3910413187914d0, -2457.9617037326980d0, -54.246284684757562d0, 500.79512237094178d0,  &
                    -1771.4461242822615d0, 3007.1993994548748d0, -2457.9617037326980d0, 776.42840231397884d0 &
                    ], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.99215674164922152d0,  &
                    1.0000000000000000d0, 0.99498743710661997d0, 0.98333216603563356d0, 0.96824583655185414d0,  &
                    0.95148591360407542d0, 0.97499604304356913d0, 0.99498743710661997d0, 1.0000000000000000d0,  &
                    0.99652172859178323d0, 0.98809481374347152d0, 0.97677102365552460d0, 0.95393920141694566d0,  &
                    0.98333216603563356d0, 0.99652172859178323d0, 1.0000000000000000d0, 0.99744571741206722d0,  &
                    0.99107124982123374d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.90905934288630952d0,  &
                    0.95148591360407542d0, 0.97677102365552460d0, 0.99107124982123374d0, 0.99804496391695696d0,  &
                    1.0000000000000000d0], shape(S))
      case(7)
        W = reshape([80.112087744414481d0, -380.83753196337295d0, 846.68851603934706d0, -1097.4660724104394d0,  &
                    854.00203553098982d0, -371.15373926820087d0, 69.379320671488060d0, -380.83753196337295d0,  &
                    2460.4826695826819d0, -6806.5629280860976d0, 10292.638409627154d0, -8930.9770566136904d0,  &
                    4194.6906070515997d0, -829.33608862620451d0, 846.68851603934706d0, -6806.5629280860976d0,  &
                    22347.744248873638d0, -38562.685975984692d0, 37023.176288616021d0, -18793.587621101400d0,  &
                    3945.6376780219634d0, -1097.4660724104394d0, 10292.638409627154d0, -38562.685975984692d0,  &
                    74347.213817597105d0, -78269.402026367752d0, 42877.575814051277d0, -9587.7268556762156d0,  &
                    854.00203553098982d0, -8930.9770566136904d0, 37023.176288616021d0, -78269.402026367752d0,  &
                    89507.134104252254d0, -52778.316794353603d0, 12594.783311461077d0, -371.15373926820087d0,  &
                    4194.6906070515997d0, -18793.587621101400d0, 42877.575814051277d0, -52778.316794353603d0,  &
                    33381.624954509396d0, -8510.6924230665845d0, 69.379320671488060d0, -829.33608862620451d0,  &
                    3945.6376780219634d0, -9587.7268556762156d0, 12594.783311461077d0, -8510.6924230665845d0,  &
                    2318.7297663139648d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0, 0.98333216603563356d0,  &
                    0.96824583655185414d0, 0.95148591360407542d0, 0.93404977361585861d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0, 0.99652172859178323d0, 0.98809481374347152d0,  &
                    0.97677102365552460d0, 0.96378881965339736d0, 0.95393920141694566d0, 0.98333216603563356d0,  &
                    0.99652172859178323d0, 1.0000000000000000d0, 0.99744571741206722d0, 0.99107124982123374d0,  &
                    0.98226460284385697d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.90905934288630952d0, 0.95148591360407542d0, 0.97677102365552460d0, 0.99107124982123374d0,  &
                    0.99804496391695696d0, 1.0000000000000000d0, 0.99845559753396829d0, 0.88712019959006128d0,  &
                    0.93404977361585861d0, 0.96378881965339736d0, 0.98226460284385697d0, 0.99305547153730200d0,  &
                    0.99845559753396829d0, 1.0000000000000000d0], shape(S))
      case(8)
        W = reshape([92.539024757823469d0, -485.62295892374391d0, 1245.3098138485498d0, -1970.1315498890267d0,  &
                    2024.2063495743746d0, -1321.8047757365766d0, 499.34651114023347d0, -83.121761910171273d0,  &
                    -485.62295892374391d0, 3539.9434708348890d0, -11511.708266971060d0, 21569.895229194808d0,  &
                    -24986.945984915284d0, 17767.671463504408d0, -7134.3810386646874d0, 1241.2366889962939d0,  &
                    1245.3098138485498d0, -11511.708266971060d0, 45076.344590889283d0, -97380.089133774367d0,  &
                    125682.09286798064d0, -97005.646137706557d0, 41462.330953470708d0, -7568.2374319869477d0,  &
                    -1970.1315498890267d0, 21569.895229194808d0, -97380.089133774367d0, 236558.01704791249d0,  &
                    -335931.86898086715d0, 280178.31947277795d0, -127521.98045884802d0, 24497.997672243826d0,  &
                    2024.2063495743746d0, -24986.945984915284d0, 125682.09286798064d0, -335931.86898086715d0,  &
                    518562.62266775075d0, -464898.16962265247d0, 225191.00178569887d0, -45642.664557060503d0,  &
                    -1321.8047757365766d0, 17767.671463504408d0, -97005.646137706557d0, 280178.31947277795d0,  &
                    -464898.16962265247d0, 445491.35949089238d0, -229343.07661604194d0, 49131.675673804595d0,  &
                    499.34651114023347d0, -7134.3810386646874d0, 41462.330953470708d0, -127521.98045884802d0,  &
                    225191.00178569887d0, -229343.07661604194d0, 125237.47296520101d0, -28390.563739949452d0,  &
                    -83.121761910171273d0, 1241.2366889962939d0, -7568.2374319869477d0, 24497.997672243826d0,  &
                    -45642.664557060503d0, 49131.675673804595d0, -28390.563739949452d0, 6814.4475643740907d0 &
                    ], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.86602540378443860d0, 0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0,  &
                    0.98333216603563356d0, 0.96824583655185414d0, 0.95148591360407542d0, 0.93404977361585861d0,  &
                    0.91651513899116799d0, 0.97499604304356913d0, 0.99498743710661997d0, 1.0000000000000000d0,  &
                    0.99652172859178323d0, 0.98809481374347152d0, 0.97677102365552460d0, 0.96378881965339736d0,  &
                    0.94991775959816649d0, 0.95393920141694566d0, 0.98333216603563356d0, 0.99652172859178323d0,  &
                    1.0000000000000000d0, 0.99744571741206722d0, 0.99107124982123374d0, 0.98226460284385697d0,  &
                    0.97192421422695907d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.98601329718326935d0, 0.90905934288630952d0, 0.95148591360407542d0, 0.97677102365552460d0,  &
                    0.99107124982123374d0, 0.99804496391695696d0, 1.0000000000000000d0, 0.99845559753396829d0,  &
                    0.99444440145743085d0, 0.88712019959006128d0, 0.93404977361585861d0, 0.96378881965339736d0,  &
                    0.98226460284385697d0, 0.99305547153730200d0, 0.99845559753396829d0, 1.0000000000000000d0,  &
                    0.99874921777190884d0, 0.86602540378443860d0, 0.91651513899116799d0, 0.94991775959816649d0,  &
                    0.97192421422695907d0, 0.98601329718326935d0, 0.99444440145743085d0, 0.99874921777190884d0,  &
                    1.0000000000000000d0], shape(S))
      case(9)
        W = reshape([102.73996023131654d0, -579.28008477226797d0, 1649.6329933089637d0, -3023.3784007460226d0,  &
                    3801.6152070879762d0, -3281.1404539213454d0, 1861.8484048267730d0, -625.49270608226868d0,  &
                    94.172577622820413d0, -579.28008477226797d0, 4626.3192102608891d0, -17024.980936140477d0,  &
                    37570.528484501476d0, -53968.531249935149d0, 51220.175635520885d0, -31105.502641416926d0,  &
                    10973.451092982890d0, -1712.0993861885026d0, 1649.6329933089637d0, -17024.980936140477d0,  &
                    76648.552042974829d0, -197491.45685094723d0, 319059.39758025052d0, -330888.18983548041d0,  &
                    214867.87633591256d0, -79759.988337705232d0, 12939.555536288681d0, -3023.3784007460226d0,  &
                    37570.528484501476d0, -197491.45685094723d0, 577420.51966726186d0, -1032806.0301229986d0,  &
                    1161968.8118061803d0, -805283.01678246027d0, 314926.85573368822d0, -53282.713938116933d0,  &
                    3801.6152070879762d0, -53968.531249935149d0, 319059.39758025052d0, -1032806.0301229986d0,  &
                    2015313.2813086673d0, -2441081.2560675377d0, 1800654.0587436187d0, -742277.93893866800d0,  &
                    131305.73650891884d0, -3281.1404539213454d0, 51220.175635520885d0, -330888.18983548041d0,  &
                    1161968.8118061803d0, -2441081.2560675377d0, 3159490.1849562805d0, -2472607.0610428583d0,  &
                    1074278.1958661408d0, -199099.57483824741d0, 1861.8484048267730d0, -31105.502641416926d0,  &
                    214867.87633591256d0, -805283.01678246027d0, 1800654.0587436187d0, -2472607.0610428583d0,  &
                    2045686.8262865648d0, -936138.89101480751d0, 182064.24577301860d0, -625.49270608226868d0,  &
                    10973.451092982890d0, -79759.988337705232d0, 314926.85573368822d0, -742277.93893866800d0,  &
                    1074278.1958661408d0, -936138.89101480751d0, 450714.89061671734d0, -92090.959139143888d0,  &
                    94.172577622820413d0, -1712.0993861885026d0, 12939.555536288681d0, -53282.713938116933d0,  &
                    131305.73650891884d0, -199099.57483824741d0, 182064.24577301860d0, -92090.959139143888d0,  &
                    19782.408284583180d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.86602540378443860d0, 0.84590516936330140d0, 0.99215674164922152d0, 1.0000000000000000d0,  &
                    0.99498743710661997d0, 0.98333216603563356d0, 0.96824583655185414d0, 0.95148591360407542d0,  &
                    0.93404977361585861d0, 0.91651513899116799d0, 0.89921841062113494d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0, 0.99652172859178323d0, 0.98809481374347152d0,  &
                    0.97677102365552460d0, 0.96378881965339736d0, 0.94991775959816649d0, 0.93564551297569798d0,  &
                    0.95393920141694566d0, 0.98333216603563356d0, 0.99652172859178323d0, 1.0000000000000000d0,  &
                    0.99744571741206722d0, 0.99107124982123374d0, 0.98226460284385697d0, 0.97192421422695907d0,  &
                    0.96064535921058791d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.98601329718326935d0, 0.97758819057930046d0, 0.90905934288630952d0, 0.95148591360407542d0,  &
                    0.97677102365552460d0, 0.99107124982123374d0, 0.99804496391695696d0, 1.0000000000000000d0,  &
                    0.99845559753396829d0, 0.99444440145743085d0, 0.98868599666425949d0, 0.88712019959006128d0,  &
                    0.93404977361585861d0, 0.96378881965339736d0, 0.98226460284385697d0, 0.99305547153730200d0,  &
                    0.99845559753396829d0, 1.0000000000000000d0, 0.99874921777190884d0, 0.99545452192223205d0,  &
                    0.86602540378443860d0, 0.91651513899116799d0, 0.94991775959816649d0, 0.97192421422695907d0,  &
                    0.98601329718326935d0, 0.99444440145743085d0, 0.99874921777190884d0, 1.0000000000000000d0,  &
                    0.99896640799254144d0, 0.84590516936330140d0, 0.89921841062113494d0, 0.93564551297569798d0,  &
                    0.96064535921058791d0, 0.97758819057930046d0, 0.98868599666425949d0, 0.99545452192223205d0,  &
                    0.99896640799254144d0, 1.0000000000000000d0], shape(S))
      case(10)
        W = reshape([110.85899098970923d0, -656.32724208730247d0, 2002.0561806262745d0, -4030.0430894678625d0,  &
                    5747.1082633813776d0, -5876.5673849153736d0, 4220.3047498562373d0, -2019.2279545723384d0,  &
                    576.83178827809900d0, -74.279463179713289d0, -656.32724208730247d0, 5580.5224445491212d0,  &
                    -22331.612043043937d0, 54995.143701985420d0, -91095.538134775474d0, 104230.69173430053d0,  &
                    -81616.918643421988d0, 41841.558030521599d0, -12662.655255751506d0, 1715.2096007504856d0,  &
                    2002.0561806262745d0, -22331.612043043937d0, 110954.23971315096d0, -323698.56673909881d0,  &
                    611814.99361075275d0, -776207.90529295243d0, 659769.38292835664d0, -361528.21759789030d0,  &
                    115656.71780639346d0, -16430.693998954230d0, -4030.0430894678625d0, 54995.143701985420d0,  &
                    -323698.56673909881d0, 1086072.8812834970d0, -2301554.0722570983d0, 3206429.9814929799d0,  &
                    -2943668.1793391723d0, 1719782.1674985890d0, -580745.53209077963d0, 86416.331611407732d0,  &
                    5747.1082633813776d0, -91095.538134775474d0, 611814.99361075275d0, -2301554.0722570983d0,  &
                    5380578.4815178504d0, -8152091.3190899873d0, 8040206.2774941670d0, -4994941.7070384044d0,  &
                    1778403.4588345746d0, -277067.37753186666d0, -5876.5673849153736d0, 104230.69173430053d0,  &
                    -776207.90529295243d0, 3206429.9814929799d0, -8152091.3190899873d0, 13309104.674707124d0,  &
                    -14025079.427185630d0, 9238846.3749270160d0, -3464500.8524018270d0, 565144.53718240885d0,  &
                    4220.3047498562373d0, -81616.918643421988d0, 659769.38292835664d0, -2943668.1793391723d0,  &
                    8040206.2774941670d0, -14025079.427185630d0, 15706401.821684649d0, -10938310.914757373d0,  &
                    4315299.1100801602d0, -737221.24287861609d0, -2019.2279545723384d0, 41841.558030521599d0,  &
                    -361528.21759789030d0, 1719782.1674985890d0, -4994941.7070384044d0, 9238846.3749270160d0,  &
                    -10938310.914757373d0, 8029013.8383061187d0, -3328353.6668580547d0, 595670.14588087215d0,  &
                    576.83178827809900d0, -12662.655255751506d0, 115656.71780639346d0, -580745.53209077963d0,  &
                    1778403.4588345746d0, -3464500.8524018270d0, 4315299.1100801602d0, -3328353.6668580547d0,  &
                    1447833.8533988285d0, -271507.14094415621d0, -74.279463179713289d0, 1715.2096007504856d0,  &
                    -16430.693998954230d0, 86416.331611407732d0, -277067.37753186666d0, 565144.53718240885d0,  &
                    -737221.24287861609d0, 595670.14588087215d0, -271507.14094415621d0, 53355.279479472047d0 &
                    ], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.86602540378443860d0, 0.84590516936330140d0, 0.82679728470768454d0, 0.99215674164922152d0,  &
                    1.0000000000000000d0, 0.99498743710661997d0, 0.98333216603563356d0, 0.96824583655185414d0,  &
                    0.95148591360407542d0, 0.93404977361585861d0, 0.91651513899116799d0, 0.89921841062113494d0,  &
                    0.88235294117647056d0, 0.97499604304356913d0, 0.99498743710661997d0, 1.0000000000000000d0,  &
                    0.99652172859178323d0, 0.98809481374347152d0, 0.97677102365552460d0, 0.96378881965339736d0,  &
                    0.94991775959816649d0, 0.93564551297569798d0, 0.92128466398761111d0, 0.95393920141694566d0,  &
                    0.98333216603563356d0, 0.99652172859178323d0, 1.0000000000000000d0, 0.99744571741206722d0,  &
                    0.99107124982123374d0, 0.98226460284385697d0, 0.97192421422695907d0, 0.96064535921058791d0,  &
                    0.94882928301683922d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.98601329718326935d0, 0.97758819057930046d0, 0.96824583655185426d0, 0.90905934288630952d0,  &
                    0.95148591360407542d0, 0.97677102365552460d0, 0.99107124982123374d0, 0.99804496391695696d0,  &
                    1.0000000000000000d0, 0.99845559753396829d0, 0.99444440145743085d0, 0.98868599666425949d0,  &
                    0.98169181562325258d0, 0.88712019959006128d0, 0.93404977361585861d0, 0.96378881965339736d0,  &
                    0.98226460284385697d0, 0.99305547153730200d0, 0.99845559753396829d0, 1.0000000000000000d0,  &
                    0.99874921777190884d0, 0.99545452192223205d0, 0.99065885080469862d0, 0.86602540378443860d0,  &
                    0.91651513899116799d0, 0.94991775959816649d0, 0.97192421422695907d0, 0.98601329718326935d0,  &
                    0.99444440145743085d0, 0.99874921777190884d0, 1.0000000000000000d0, 0.99896640799254144d0,  &
                    0.99621210759909562d0, 0.84590516936330140d0, 0.89921841062113494d0, 0.93564551297569798d0,  &
                    0.96064535921058791d0, 0.97758819057930046d0, 0.98868599666425949d0, 0.99545452192223205d0,  &
                    0.99896640799254144d0, 1.0000000000000000d0, 0.99913156735681652d0, 0.82679728470768454d0,  &
                    0.88235294117647056d0, 0.92128466398761111d0, 0.94882928301683922d0, 0.96824583655185426d0,  &
                    0.98169181562325258d0, 0.99065885080469862d0, 0.99621210759909562d0, 0.99913156735681652d0,  &
                    1.0000000000000000d0], shape(S))
      case(11)
        W = reshape([117.21235214772904d0, -714.35806297792158d0, 2245.1643321280740d0, -4610.3933792100061d0,  &
                    6511.4009202362895d0, -6150.8302224825920d0, 3391.5937860140148d0, -436.01742880721537d0,  &
                    -738.92275647372799d0, 484.31715587976799d0, -98.454105559226775d0, -714.35806297792158d0,  &
                    6269.8952732304706d0, -25875.069939432633d0, 65110.259187545853d0, -107732.04536419440d0,  &
                    117390.86907425059d0, -78379.654709546943d0, 23457.309598262604d0, 5607.1440832265725d0,  &
                    -6684.8214165052414d0, 1550.5418761193073d0, 2245.1643321280740d0, -25875.069939432633d0,  &
                    132599.60954916809d0, -396361.31134008802d0, 756826.91023577540d0, -944113.58367922041d0,  &
                    746475.35184183146d0, -329819.32141719793d0, 40631.572575033351d0, 26062.477344374791d0,  &
                    -8671.4089555561513d0, -4610.3933792100061d0, 65110.259187545853d0, -396361.31134008802d0,  &
                    1372463.1491815776d0, -2988502.0270795608d0, 4244631.4027743703d0, -3916425.4633714431d0,  &
                    2232511.7325943527d0, -677804.31608534034d0, 55391.291470506068d0, 13595.784351342649d0,  &
                    6511.4009202362895d0, -107732.04536419440d0, 756826.91023577540d0, -2988502.0270795608d0,  &
                    7374201.5953415921d0, -11897330.128776217d0, 12685063.180279443d0, -8765094.4377386160d0,  &
                    3697769.2054832382d0, -829684.91928406083d0, 67971.552751104493d0, -6150.8302224825920d0,  &
                    117390.86907425059d0, -944113.58367922041d0, 4244631.4027743703d0, -11897330.128776217d0,  &
                    21851270.467823889d0, -26754279.226484302d0, 21633158.682814363d0, -11077651.903835505d0,  &
                    3246049.5859478437d0, -412975.14465735806d0, 3391.5937860140148d0, -78379.654709546943d0,  &
                    746475.35184183146d0, -3916425.4633714431d0, 12685063.180279443d0, -26754279.226484302d0,  &
                    37494848.639549591d0, -34689165.786953673d0, 20377863.358143974d0, -6891923.6598007092d0,  &
                    1022531.8799494690d0, -436.01742880721537d0, 23457.309598262604d0, -329819.32141719793d0,  &
                    2232511.7325943527d0, -8765094.4377386160d0, 21633158.682814363d0, -34689165.786953673d0,  &
                    36151306.746357322d0, -23647318.572119914d0, 8826384.4847903624d0, -1434984.5992588799d0,  &
                    -738.92275647372799d0, 5607.1440832265725d0, 40631.572575033351d0, -677804.31608534034d0,  &
                    3697769.2054832382d0, -11077651.903835505d0, 20377863.358143974d0, -23647318.572119914d0,  &
                    16913806.476479985d0, -6818990.6947852084d0, 1186826.9821223991d0, 484.31715587976799d0,  &
                    -6684.8214165052414d0, 26062.477344374791d0, 55391.291470506068d0, -829684.91928406083d0,  &
                    3246049.5859478437d0, -6891923.6598007092d0, 8826384.4847903624d0, -6818990.6947852084d0,  &
                    2933523.2297482090d0, -540611.16396153928d0, -98.454105559226775d0, 1550.5418761193073d0,  &
                    -8671.4089555561513d0, 13595.784351342649d0, 67971.552751104493d0, -412975.14465735806d0,  &
                    1022531.8799494690d0, -1434984.5992588799d0, 1186826.9821223991d0, -540611.16396153928d0,  &
                    104864.79597292397d0], shape(W))
        S = reshape([1.0000000000000000d0, 0.99215674164922152d0, 0.97499604304356913d0,  &
                    0.95393920141694566d0, 0.93154097872359987d0, 0.90905934288630952d0, 0.88712019959006128d0,  &
                    0.86602540378443860d0, 0.84590516936330140d0, 0.82679728470768454d0, 0.80868982852161886d0,  &
                    0.99215674164922152d0, 1.0000000000000000d0, 0.99498743710661997d0, 0.98333216603563356d0,  &
                    0.96824583655185414d0, 0.95148591360407542d0, 0.93404977361585861d0, 0.91651513899116799d0,  &
                    0.89921841062113494d0, 0.88235294117647056d0, 0.86602540378443871d0, 0.97499604304356913d0,  &
                    0.99498743710661997d0, 1.0000000000000000d0, 0.99652172859178323d0, 0.98809481374347152d0,  &
                    0.97677102365552460d0, 0.96378881965339736d0, 0.94991775959816649d0, 0.93564551297569798d0,  &
                    0.92128466398761111d0, 0.90703620734810975d0, 0.95393920141694566d0, 0.98333216603563356d0,  &
                    0.99652172859178323d0, 1.0000000000000000d0, 0.99744571741206722d0, 0.99107124982123374d0,  &
                    0.98226460284385697d0, 0.97192421422695907d0, 0.96064535921058791d0, 0.94882928301683922d0,  &
                    0.93674969975975964d0, 0.93154097872359987d0, 0.96824583655185414d0, 0.98809481374347152d0,  &
                    0.99744571741206722d0, 1.0000000000000000d0, 0.99804496391695696d0, 0.99305547153730200d0,  &
                    0.98601329718326935d0, 0.97758819057930046d0, 0.96824583655185426d0, 0.95831484749990992d0,  &
                    0.90905934288630952d0, 0.95148591360407542d0, 0.97677102365552460d0, 0.99107124982123374d0,  &
                    0.99804496391695696d0, 1.0000000000000000d0, 0.99845559753396829d0, 0.99444440145743085d0,  &
                    0.98868599666425949d0, 0.98169181562325258d0, 0.97383114934675230d0, 0.88712019959006128d0,  &
                    0.93404977361585861d0, 0.96378881965339736d0, 0.98226460284385697d0, 0.99305547153730200d0,  &
                    0.99845559753396829d0, 1.0000000000000000d0, 0.99874921777190884d0, 0.99545452192223205d0,  &
                    0.99065885080469862d0, 0.98476101329618471d0, 0.86602540378443860d0, 0.91651513899116799d0,  &
                    0.94991775959816649d0, 0.97192421422695907d0, 0.98601329718326935d0, 0.99444440145743085d0,  &
                    0.99874921777190884d0, 1.0000000000000000d0, 0.99896640799254144d0, 0.99621210759909562d0,  &
                    0.99215674164922152d0, 0.84590516936330140d0, 0.89921841062113494d0, 0.93564551297569798d0,  &
                    0.96064535921058791d0, 0.97758819057930046d0, 0.98868599666425949d0, 0.99545452192223205d0,  &
                    0.99896640799254144d0, 1.0000000000000000d0, 0.99913156735681652d0, 0.99679486355016889d0,  &
                    0.82679728470768454d0, 0.88235294117647056d0, 0.92128466398761111d0, 0.94882928301683922d0,  &
                    0.96824583655185426d0, 0.98169181562325258d0, 0.99065885080469862d0, 0.99621210759909562d0,  &
                    0.99913156735681652d0, 1.0000000000000000d0, 0.99926008128973698d0, 0.80868982852161886d0,  &
                    0.86602540378443871d0, 0.90703620734810975d0, 0.93674969975975964d0, 0.95831484749990992d0,  &
                    0.97383114934675230d0, 0.98476101329618471d0, 0.99215674164922152d0, 0.99679486355016889d0,  &
                    0.99926008128973698d0, 1.0000000000000000d0], shape(S))
      case(12:)
        write(*,*) "Bad value of alpha_max"
        stop
    end select

  return
  end subroutine
!**************************************************************************


!-----------------------------------------------------------------------------------
! OPERATOR FORM OF RADIAL EXP. COEFF. STARTS HERE
!-----------------------------------------------------------------------------------
!
!**************************************************************************
!
! This subroutine returns the radial expansion coefficients using the
! polynomial basis set and a polynomial piecewise representation for 
! the atomic density
!
  subroutine get_radial_expansion_coefficients_poly3operator(n_sites, n_neigh, rjs_in, alpha_max, &
                                                             rcut_soft_in, rcut_hard_in, atom_sigma_in, &
                                                             atom_sigma_scaling, amplitude_scaling, W, &
                                                             scaling_mode, mask, radial_enhancement, &
                                                             do_derivatives, do_central, central_weight, &
                                                             exp_coeff, exp_coeff_der)

    implicit none

    integer, intent(in) :: alpha_max, n_neigh(:), n_sites, radial_enhancement
    real(dp), intent(in) :: rcut_soft_in, rcut_hard_in, rjs_in(:), atom_sigma_in, atom_sigma_scaling
    real(dp), intent(in) :: amplitude_scaling, central_weight
    real(dp) :: rcut_soft, rcut_hard, atom_sigma, atom_sigma_scaled, amplitude
    logical, intent(in) :: mask(:), do_derivatives, do_central
    character(*), intent(in) :: scaling_mode
!
    integer :: i, j, k
    real(dp) :: pi, rj, s2, atom_width, atom_width_scaling, filter_width, x
    real(dp), allocatable :: A(:,:)
    real(dp) :: W(:,:)
!   derivatives
    real(dp) :: amplitude_der
!   Results will be stored in exp_coeff, which is an array of dimension (alpha_max, n_atom_pairs)
    real(dp) :: exp_coeff(:,:), exp_coeff_der(:,:)
    logical, save :: print_basis = .false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(dp), allocatable :: rjs(:), lim_soft_array(:,:), atom_sigma_scaleds(:), s2s(:), atom_widths(:), &
                           I0_array(:,:,:), g_aux_left_array(:,:,:), g_aux_right_array(:,:,:), &
                           M_left_array(:,:,:), M_right_array(:,:,:), exp_coeff_soft_array(:,:), &
                           amplitudes(:), I_left_array(:,:), I_right_array(:,:), B(:,:), &
                           exp_coeff_buffer_array(:,:), lim_buffer_array(:,:)
!!!! derivatives
    real(dp), allocatable :: g_aux_left_der_array(:,:,:), g_aux_right_der_array(:,:,:), &
                           M_left_der_array(:,:,:), M_right_der_array(:,:,:), B_der(:,:), &
                           exp_coeff_soft_der_array(:,:), amplitudes_der(:), &
                           I_left_der_array(:,:), I_right_der_array(:,:), exp_coeff_buffer_der_array(:,:)

    real(dp) :: vect(1:4)
    integer, allocatable :: rjs_idx(:)
!    integer, allocatable ::  soft_weights(:), buffer_weights(:)
    integer :: nn, k2
    real(dp), save :: elapsed_time = 0.d0
    real(dp) :: t1, t2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   This is for debugging. It prints the basis set to plot it with Gnuplot (gfortran only)
    if( print_basis )then
      print_basis = .false.
      write(*,*) "p(x,n,rcut) = (1.-x/rcut)**(n+2) / sqrt( rcut / (2.*n+5.) ) "
      do j=1, alpha_max
        write(*,"(A,I0,A)",advance="no") "p", j, "(x) = "
        do i = 1, alpha_max
          write(*,"(A,I2,A,F16.10,A,E16.8,A)",advance="no") "p(x,", i, ",", rcut_hard_in, ") *", W(j,i), "+"
        end do
        write(*,*) "0."
      end do
    end if

!   *************** New basis *******************
!   Normalized with rcut_hard to avoid instabilities
    rcut_soft = rcut_soft_in/rcut_hard_in
    rcut_hard = 1.d0
    atom_sigma = atom_sigma_in/rcut_hard_in
!   *********************************************
    pi = dacos(-1.d0)
    filter_width = 2.d0*sqrt(2.d0*log(2.d0))*(rcut_hard - rcut_soft)
    exp_coeff = 0.d0
    if( do_derivatives )then
      exp_coeff_der = 0.d0
    end if

    allocate( A(1:alpha_max, 1:7) )
    call get_constant_poly_coeff(alpha_max, rcut_hard, A)

! NOTE: the vectorization can potentially be done for all atoms at one, where nn = size(rjs_in), then only
! obvious thing to be careful about I can think of is handling of the central atoms in the flattened array
    k = 0
    do i = 1, n_sites
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   check number of atoms within soap cutoff
!     nn = count( rjs_in(k+1:k+n_neigh(i)) < rcut_hard_in .and. mask(k+1:k+n_neigh(i)) )
!     if( nn < 1 )then
!       cycle
!     end if


!!! SOFT REGION
!     Temporarily allocate this as oversized array:
      allocate( atom_widths(1:n_neigh(i)) )
      atom_widths = 2.d0*sqrt(2.d0*log(2.d0))*(atom_sigma_in + atom_sigma_scaling*rjs_in(k+1:k+n_neigh(i)))

!     count number of atoms fully within soft region
      nn = count( rjs_in(k+1:k+n_neigh(i)) - atom_widths < rcut_soft_in .and. mask(k+1:k+n_neigh(i)) )
      if( nn > 0 )then
!       These need to be allocated immediately
        allocate( rjs(1:nn) )
        allocate( rjs_idx(1:nn) )

!       make assignment between original neighbors and reduced neighbors
        k2 = 1
        do j = 1, n_neigh(i)
!         if( rjs_in(k+j) < rcut_hard_in .and. mask(k+j) )then
          if( rjs_in(k+j) - atom_widths(j) < rcut_soft_in .and. mask(k+j) )then
            rjs(k2) = rjs_in(k+j)
            rjs_idx(k2) = k+j
            k2 = k2 + 1
          end if
        end do

        deallocate( atom_widths )

!       allocate arrays that depend on the number of neighbors inside soap cutoff
        allocate( atom_sigma_scaleds(1:nn) )
        allocate( s2s(1:nn) )
        allocate( atom_widths(1:nn) )
        allocate( lim_soft_array(1:nn, 1:3) )
!        allocate( soft_weights(1:nn) )
!        allocate( buffer_weights(1:nn) )
!       we might be able to get some extra speed by swapping the order of some of these dimensions
        allocate( I0_array(1:nn, 1:alpha_max + 4, 1:3) )
        allocate( g_aux_left_array(1:nn, 1:alpha_max, 1:2) )
        allocate( g_aux_right_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
        allocate( M_left_array(1:nn, 1:4, 1:2) )
        allocate( M_right_array(1:nn, 1:4, 2:3) ) ! note the 2:3 bounds here
        allocate( I_left_array(1:nn, 1:alpha_max) )
        allocate( I_right_array(1:nn, 1:alpha_max) )
        allocate( exp_coeff_soft_array(1:nn, 1:alpha_max) )
!        allocate( exp_coeff_buffer_array(1:nn, 1:alpha_max) )
        allocate( amplitudes(1:nn) )
        allocate( amplitudes_der(1:nn) )
        allocate( B(1:7, 1:nn) )

        if( do_derivatives )then
          allocate( g_aux_left_der_array(1:nn, 1:alpha_max, 1:2) )
          allocate( g_aux_right_der_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
          allocate( M_left_der_array(1:nn, 1:4, 1:2) )
          allocate( M_right_der_array(1:nn, 1:4, 2:3) ) ! note the 2:3 bounds here
          allocate( I_left_der_array(1:nn, 1:alpha_max) )
          allocate( I_right_der_array(1:nn, 1:alpha_max) )
          allocate( exp_coeff_soft_der_array(1:nn, 1:alpha_max) )
        end if

!       define atom parameters in vector form
!   **************** New basis ******************
        rjs = rjs/rcut_hard_in
!   *********************************************
        atom_sigma_scaleds = atom_sigma + atom_sigma_scaling*rjs
        s2s = atom_sigma_scaleds**2
        atom_widths = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaleds
        atom_width_scaling = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaling
!   ~~~~~~~~~~~~~~~ Amplitude ~~~~~~~~~~~~~~~~~~~~~~~
        if( scaling_mode == "polynomial" )then
!         WARNING: the 1/atom_sigma_angular^2 term is missing from these amplitudes and needs to
!                  be taken into account in the corresponding part of the code.
!       
!         WARNING2: These expressions here already assume rcut_hard = 1., so this parameter is missing
!                   from the expressions
          if( amplitude_scaling == 0.d0 )then
            amplitudes = 1.d0 / atom_sigma_scaleds
            amplitudes_der = - atom_sigma_scaling / s2s 
!            else if( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 <= 1.d-10 )then ! FIX THESE !!!!!!!!!!!!!!!!!!!!!!!!!!
!              amplitudes = 0.d0
!              amplitude_ders = 0.d0
          else
            if( amplitude_scaling == 1.d0 )then
              amplitudes = 1.d0 / atom_sigma_scaleds * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )
              amplitudes_der = 6.d0 / atom_sigma_scaleds * (rjs**2 - rjs) &
                                - atom_sigma_scaling / atom_sigma_scaleds * amplitudes
            else
              amplitudes = 1.d0 / atom_sigma_scaleds * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )**amplitude_scaling
              amplitudes_der = 6.d0*amplitude_scaling / atom_sigma_scaleds * (rjs**2 - rjs) &
                                * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )**(amplitude_scaling - 1.d0) &
                                - atom_sigma_scaling / atom_sigma_scaleds * amplitudes
            end if
          end if
        end if
!       The central atom needs to be scaled by central_weight
        if( size(amplitudes) > 0 )then
          if( rjs_idx(1) == k+1 )then
            amplitudes(1) = central_weight * amplitudes(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            amplitudes_der(1) = central_weight * amplitudes_der(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          end if
        end if

!       The radial enhancement adds a scaling corresponding to the integral of a Gaussian at the position
!       of atom j.
        if( radial_enhancement == 1 )then
          amplitudes_der = amplitudes * ( 1.d0 + dsqrt(2.d0/pi)*atom_sigma_scaling ) + &
                           amplitudes_der * ( rjs + dsqrt(2.d0/pi)*atom_sigma_scaleds )
          amplitudes = amplitudes * ( rjs + dsqrt(2.d0/pi)*atom_sigma_scaleds )
        else if( radial_enhancement == 2 )then
          amplitudes_der = amplitudes*( 2.d0*rjs + 2.d0*atom_sigma_scaleds*atom_sigma_scaling + &
                           dsqrt(8.d0/pi)*atom_sigma_scaleds + dsqrt(8.d0/pi)*rjs*atom_sigma_scaling ) + &
                           amplitudes_der*( rjs**2 + s2s + dsqrt(8.d0/pi)*atom_sigma_scaleds*rjs )
          amplitudes = amplitudes * ( rjs**2 + s2s + dsqrt(8.d0/pi)*atom_sigma_scaleds*rjs )
        end if     
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!       NOTE: this below is design to be used with arrays that contain all the atoms. In some cases, the 
!       left_weights or right_weights might contain lots of zeros, in which case a speed up of close to x2
!       could be achieved by not doing all atoms but only those for which the weights are 1. 
!       One could refactor this code with the nn construction to have, e.g., nn_soft and nn_buffer
!       define arrays that give contribution of atoms within left and right of rcut_soft
!       Note that these are integer arrays (0|1 = .false.|.true.)
!        soft_weights = rjs - atom_widths < rcut_soft
!        nn_soft = count(soft_weights)
!        buffer_weights = rcut_soft < rcut_hard .and. rjs + atom_widths > rcut_soft
!        nn_buffer = count(buffer_weights)

!       handle the central atom contribution; 
!       NOTE: I'M NOT SURE ABOUT THE NEED FOR THIS DO_CENTRAL TAG HERE, IT SHOULD BE HANDLED VIA THE CENTRAL_WEIGHT
        if( .not. do_central )then
          if( size(amplitudes) > 0 )then
            if( rjs_idx(1) == k+1 )then
              amplitudes(1) = 0.d0
            end if
          end if
        end if

!       integration limits inside r_soft
        lim_soft_array(:, 1) = max( 0.d0, rjs - atom_widths )      ! lower limit left
        lim_soft_array(:, 2) = min( rjs, rcut_soft )               ! upper limit left / lower limit right
        lim_soft_array(:, 3) = min( rcut_soft, rjs + atom_widths ) ! upper limit right

!       dimensions of M_radial_poly_array are (1:nn, 1:alpha_max + 4, 1:3)
        I0_array = M_radial_poly_array(lim_soft_array, alpha_max + 4, rcut_hard)
        g_aux_left_array(1:nn, 1:4, 1:2) = g_aux_array(lim_soft_array(:, 1:2), rjs(:), atom_widths(:), "left")
        g_aux_right_array(1:nn, 1:4, 2:3) = g_aux_array(lim_soft_array(:, 2:3), rjs(:), atom_widths(:), "right")

        vect = [-1.d0, -1.d0, -2.d0, -6.d0]
        do j = 1, 4
          M_left_array(1:nn, j, 1) = I0_array(1:nn, j, 1) * g_aux_left_array(1:nn, j, 1) * vect(j)
        end do
        do j = 1, 4
          M_left_array(1:nn, j, 2) = I0_array(1:nn, j, 2) * g_aux_left_array(1:nn, j, 2) * vect(j)
        end do
        do j = 1, 4
          M_right_array(1:nn, j, 2) = I0_array(1:nn, j, 2) * g_aux_right_array(1:nn, j, 2) * vect(j)
        end do
        do j = 1, 4
          M_right_array(1:nn, j, 3) = I0_array(1:nn, j, 3) * g_aux_right_array(1:nn, j, 3) * vect(j)
        end do

        I_left_array(1:nn, 1:alpha_max) = matmul( M_left_array(1:nn, 1:4, 2), transpose(A(1:alpha_max, 1:4)) ) * &
                                          I0_array(1:nn, 5:alpha_max + 4, 2) - &
                                          matmul( M_left_array(1:nn, 1:4, 1), transpose(A(1:alpha_max, 1:4)) ) * &
                                          I0_array(1:nn, 5:alpha_max + 4, 1)
        I_right_array(1:nn, 1:alpha_max) = matmul( M_right_array(1:nn, 1:4, 3), transpose(A(1:alpha_max, 1:4)) ) * &
                                           I0_array(1:nn, 5:alpha_max + 4, 3) - &
                                           matmul( M_right_array(1:nn, 1:4, 2), transpose(A(1:alpha_max, 1:4)) ) * &
                                           I0_array(1:nn, 5:alpha_max + 4, 2)
        exp_coeff_soft_array = I_left_array + I_right_array

        do j = 1, nn
          k2 = rjs_idx(j)
          exp_coeff(1:alpha_max, k2) = amplitudes(j) * exp_coeff_soft_array(j, 1:alpha_max)
        end do

!       Contributions to derivatives
        if( do_derivatives )then
          g_aux_left_der_array(1:nn, 1:4, 1:2) = g_aux_der_array(lim_soft_array(:, 1:2), rjs(:), atom_widths(:), &
                                                                 atom_width_scaling, "left")
          g_aux_right_der_array(1:nn, 1:4, 2:3) = g_aux_der_array(lim_soft_array(:, 2:3), rjs(:), atom_widths(:), &
                                                                  atom_width_scaling, "right")
          do j = 1, 4
            M_left_der_array(1:nn, j, 1) = I0_array(1:nn, j, 1) * g_aux_left_der_array(1:nn, j, 1) * vect(j)
          end do
          do j = 1, 4
            M_left_der_array(1:nn, j, 2) = I0_array(1:nn, j, 2) * g_aux_left_der_array(1:nn, j, 2) * vect(j)
          end do
          do j = 1, 4
            M_right_der_array(1:nn, j, 2) = I0_array(1:nn, j, 2) * g_aux_right_der_array(1:nn, j, 2) * vect(j)
          end do
          do j = 1, 4
            M_right_der_array(1:nn, j, 3) = I0_array(1:nn, j, 3) * g_aux_right_der_array(1:nn, j, 3) * vect(j)
          end do

          I_left_der_array(1:nn, 1:alpha_max) = matmul( M_left_der_array(1:nn, 1:4, 2), transpose(A(1:alpha_max, 1:4)) ) * &
                                                I0_array(1:nn, 5:alpha_max + 4, 2) - &
                                                matmul( M_left_der_array(1:nn, 1:4, 1), transpose(A(1:alpha_max, 1:4)) ) * &
                                                I0_array(1:nn, 5:alpha_max + 4, 1)
          I_right_der_array(1:nn, 1:alpha_max) = matmul( M_right_der_array(1:nn, 1:4, 3), transpose(A(1:alpha_max, 1:4)) ) * &
                                                 I0_array(1:nn, 5:alpha_max + 4, 3) - &
                                                 matmul( M_right_der_array(1:nn, 1:4, 2), transpose(A(1:alpha_max, 1:4)) ) * &
                                                 I0_array(1:nn, 5:alpha_max + 4, 2)
          exp_coeff_soft_der_array = I_left_der_array + I_right_der_array

          do j = 1, nn
            k2 = rjs_idx(j)
            exp_coeff_der(1:alpha_max, k2) = amplitudes(j) * exp_coeff_soft_der_array(j, 1:alpha_max) + &
                                             amplitudes_der(j) * exp_coeff_soft_array(j, 1:alpha_max)
          end do

          deallocate( g_aux_left_der_array, g_aux_right_der_array, M_left_der_array, M_right_der_array, & 
                      I_left_der_array, I_right_der_array, exp_coeff_soft_der_array )

        end if

!        deallocate( soft_weights, buffer_weights )
        deallocate( rjs, rjs_idx, atom_sigma_scaleds, s2s, atom_widths, lim_soft_array, &
                    I0_array, g_aux_left_array, g_aux_right_array, M_left_array, M_right_array, &
                    I_left_array, I_right_array, exp_coeff_soft_array, amplitudes, amplitudes_der, B)
      else
        deallocate( atom_widths )
      end if


!!! BUFFER REGION
!     Temporarily allocate this as oversized array:
      allocate( atom_widths(1:n_neigh(i)) )
      atom_widths = 2.d0*sqrt(2.d0*log(2.d0))*(atom_sigma_in + atom_sigma_scaling*rjs_in(k+1:k+n_neigh(i)))

!     count number of atoms within buffer region
      nn = count( rcut_soft_in < rcut_hard_in .and. mask(k+1:k+n_neigh(i)) .and. &
                  rjs_in(k+1:k+n_neigh(i)) + atom_widths > rcut_soft_in )
      if( nn > 0 )then
!       These need to be allocated immediately
        allocate( rjs(1:nn) )
        allocate( rjs_idx(1:nn) )
!       make assignment between original neighbors and reduced neighbors
        k2 = 1
        do j = 1, n_neigh(i)
!          if( rjs_in(k+j) < rcut_hard_in .and. mask(k+j) )then
          if( rcut_soft_in < rcut_hard_in .and. rjs_in(k+j) + atom_widths(j) > rcut_soft_in .and. mask(k+j) )then
            rjs(k2) = rjs_in(k+j)
            rjs_idx(k2) = k+j
            k2 = k2 + 1
          end if
        end do

        deallocate( atom_widths )

!       allocate arrays that depend on the number of neighbors inside soap cutoff
        allocate( atom_sigma_scaleds(1:nn) )
        allocate( s2s(1:nn) )
        allocate( atom_widths(1:nn) )
        allocate( lim_buffer_array(1:nn, 1:3) )
!        allocate( soft_weights(1:nn) )
!        allocate( buffer_weights(1:nn) )
!       we might be able to get some extra speed by swapping the order of some of these dimensions
        allocate( I0_array(1:nn, 1:max(7, alpha_max + 4), 1:3) )
        allocate( g_aux_left_array(1:nn, 1:alpha_max, 1:2) )
        allocate( g_aux_right_array(1:nn, 1:alpha_max, 2:3) ) ! note the 2:3 bounds here
        allocate( M_left_array(1:nn, 1:7, 1:2) )
        allocate( M_right_array(1:nn, 1:7, 2:3) ) ! note the 2:3 bounds here
        allocate( I_left_array(1:nn, 1:alpha_max) )
        allocate( I_right_array(1:nn, 1:alpha_max) )
!        allocate( exp_coeff_soft_array(1:nn, 1:alpha_max) )
        allocate( exp_coeff_buffer_array(1:nn, 1:alpha_max) )
        allocate( amplitudes(1:nn) )
        allocate( amplitudes_der(1:nn) )
        allocate( B(1:7, 1:nn) )

        if( do_derivatives )then
          allocate( B_der(1:7, 1:nn) )
          allocate( M_left_der_array(1:nn, 1:7, 1:2) )
          allocate( M_right_der_array(1:nn, 1:7, 2:3) ) ! note the 2:3 bounds here
          allocate( I_left_der_array(1:nn, 1:alpha_max) )
          allocate( I_right_der_array(1:nn, 1:alpha_max) )
          allocate( exp_coeff_buffer_der_array(1:nn, 1:alpha_max) )
        end if

!       define atom parameters in vector form
!   **************** New basis ******************
        rjs = rjs/rcut_hard_in
!   *********************************************
        atom_sigma_scaleds = atom_sigma + atom_sigma_scaling*rjs
        s2s = atom_sigma_scaleds**2
        atom_widths = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaleds
        atom_width_scaling = 2.d0*sqrt(2.d0*log(2.d0))*atom_sigma_scaling
!   ~~~~~~~~~~~~~~~ Amplitude ~~~~~~~~~~~~~~~~~~~~~~~
        if( scaling_mode == "polynomial" )then
!         WARNING: the 1/atom_sigma_angular^2 term is missing from these amplitudes and needs to
!                  be taken into account in the corresponding part of the code.
!         WARNING2: These expressions here already assume rcut_hard = 1., so this parameter is missing
!                   from the expressions
          if( amplitude_scaling == 0.d0 )then
            amplitudes = 1.d0 / atom_sigma_scaleds
            amplitudes_der = - atom_sigma_scaling / s2s 
!         else if( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 <= 1.d-10 )then ! FIX THESE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           amplitudes = 0.d0
!           amplitude_ders = 0.d0
          else
            if( amplitude_scaling == 1.d0 )then
              amplitudes = 1.d0 / atom_sigma_scaleds * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )
              amplitudes_der = 6.d0 / atom_sigma_scaleds * (rjs**2 - rjs) &
                              - atom_sigma_scaling / atom_sigma_scaleds * amplitudes
            else
              amplitudes = 1.d0 / atom_sigma_scaleds * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )**amplitude_scaling
              amplitudes_der = 6.d0*amplitude_scaling / atom_sigma_scaleds * (rjs**2 - rjs) &
                               * ( 1.d0 + 2.d0*rjs**3 - 3.d0*rjs**2 )**(amplitude_scaling - 1.d0) &
                               - atom_sigma_scaling / atom_sigma_scaleds * amplitudes
            end if
          end if
        end if
!       The central atom needs to be scaled by central_weight
        if( size(amplitudes) > 0 )then
          if( rjs_idx(1) == k+1 )then
            amplitudes(1) = central_weight * amplitudes(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            amplitudes_der(1) = central_weight * amplitudes_der(1) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          end if
        end if
!       The radial enhancement adds a scaling corresponding to the integral of a Gaussian at the position
!       of atom j.
        if( radial_enhancement == 1 )then
          amplitudes_der = amplitudes * ( 1.d0 + dsqrt(2.d0/pi)*atom_sigma_scaling ) + &
                           amplitudes_der * ( rjs + dsqrt(2.d0/pi)*atom_sigma_scaleds )
          amplitudes = amplitudes * ( rjs + dsqrt(2.d0/pi)*atom_sigma_scaleds )
        else if( radial_enhancement == 2 )then
          amplitudes_der = amplitudes*( 2.d0*rjs + 2.d0*atom_sigma_scaleds*atom_sigma_scaling + &
                              dsqrt(8.d0/pi)*atom_sigma_scaleds + dsqrt(8.d0/pi)*rjs*atom_sigma_scaling ) + &
                           amplitudes_der*( rjs**2 + s2s + dsqrt(8.d0/pi)*atom_sigma_scaleds*rjs )
          amplitudes = amplitudes * ( rjs**2 + s2s + dsqrt(8.d0/pi)*atom_sigma_scaleds*rjs )
        end if     
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!       NOTE: this below is design to be used with arrays that contain all the atoms. In some cases, the 
!       left_weights or right_weights might contain lots of zeros, in which case a speed up of close to x2
!       could be achieved by not doing all atoms but only those for which the weights are 1. 
!       One could refactor this code with the nn construction to have, e.g., nn_soft and nn_buffer
!       define arrays that give contribution of atoms within left and right of rcut_soft
!       Note that these are integer arrays (0|1 = .false.|.true.)
!        soft_weights = rjs - atom_widths < rcut_soft
!        nn_soft = count(soft_weights)
!        buffer_weights = rcut_soft < rcut_hard .and. rjs + atom_widths > rcut_soft
!        nn_buffer = count(buffer_weights)

!       handle the central atom contribution; 
!       NOTE: I'M NOT SURE ABOUT THE NEED FOR THIS DO_CENTRAL TAG HERE, IT SHOULD BE HANDLED VIA THE CENTRAL_WEIGHT
        if( .not. do_central )then
          if( size(amplitudes) > 0 )then
            if( rjs_idx(1) == k+1 )then
              amplitudes(1) = 0.d0
            end if
          end if
        end if

!       Atoms within the buffer region
        lim_buffer_array(:, 1) = max( rcut_soft, rjs - atom_widths ) ! lower limit left
        lim_buffer_array(:, 2) = max( rjs, rcut_soft )               ! upper limit left / lower limit right
        lim_buffer_array(:, 3) = min( rcut_hard, rjs + atom_widths ) ! upper limit right

        I0_array = M_radial_poly_array(lim_buffer_array, max(7, alpha_max + 4), rcut_hard)

        call get_constant_poly_filter_coeff_array(rjs, atom_widths, rcut_soft, filter_width, 'left', B)

!       We should try to figure out a more "vectorized way" of doing this
        do k2 = 1, nn
          M_left_array(k2, 1:7, 1) = matmul( -B(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 1), 6) ) * &
                                     I0_array(k2, 1:7, 1)
          M_left_array(k2, 1:7, 2) = matmul( -B(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 2), 6) ) * &
                                     I0_array(k2, 1:7, 2)
        end do

        I_left_array(1:nn, 1:alpha_max) = matmul( M_left_array(1:nn, 1:7, 2), transpose(A(1:alpha_max, 1:7)) ) * &
                                          I0_array(1:nn, 5:alpha_max + 4, 2) - &
                                          matmul( M_left_array(1:nn, 1:7, 1), transpose(A(1:alpha_max, 1:7)) ) * &
                                          I0_array(1:nn, 5:alpha_max + 4, 1)

        call get_constant_poly_filter_coeff_array(rjs, atom_widths, rcut_soft, filter_width, 'right', B)

!       We should try to figure out a more "vectorized way" of doing this
        do k2 = 1, nn
          M_right_array(k2, 1:7, 2) = matmul( -B(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 2), 6) ) * &
                                      I0_array(k2, 1:7, 2)
          M_right_array(k2, 1:7, 3) = matmul( -B(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 3), 6) ) * &
                                      I0_array(k2, 1:7, 3)
        end do

        I_right_array(1:nn, 1:alpha_max) = matmul( M_right_array(1:nn, 1:7, 3), transpose(A(1:alpha_max, 1:7)) ) * &
                                           I0_array(1:nn, 5:alpha_max + 4, 3) - &
                                           matmul( M_right_array(1:nn, 1:7, 2), transpose(A(1:alpha_max, 1:7)) ) * &
                                           I0_array(1:nn, 5:alpha_max + 4, 2)
        exp_coeff_buffer_array = I_left_array + I_right_array

        do j = 1, nn
          k2 = rjs_idx(j)
          exp_coeff(1:alpha_max, k2) = exp_coeff(1:alpha_max, k2) + &
                                       amplitudes(j) * exp_coeff_buffer_array(j, 1:alpha_max)
        end do

!       Contribution to derivatives
        if( do_derivatives )then
          call get_constant_poly_filter_coeff_der_array(rjs, atom_widths, atom_width_scaling, rcut_soft, &
                                                        filter_width, 'left', B_der)
!         We should try to figure out a more "vectorized way" of doing this
          do k2 = 1, nn
            M_left_der_array(k2, 1:7, 1) = I0_array(k2, 1:7, 1) * &
                                           matmul( -B_der(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 1), 6) )
            M_left_der_array(k2, 1:7, 2) = I0_array(k2, 1:7, 2) * &
                                           matmul( -B_der(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 2), 6) )
          end do
          I_left_der_array(1:nn, 1:alpha_max) = matmul( M_left_der_array(1:nn, 1:7, 2), transpose(A(1:alpha_max, 1:7)) ) * &
                                                I0_array(1:nn, 5:alpha_max + 4, 2) - &
                                                matmul( M_left_der_array(1:nn, 1:7, 1), transpose(A(1:alpha_max, 1:7)) ) * &
                                                I0_array(1:nn, 5:alpha_max + 4, 1)

          call get_constant_poly_filter_coeff_der_array(rjs, atom_widths, atom_width_scaling, rcut_soft, &
                                                        filter_width, 'right', B_der)
!         We should try to figure out a more "vectorized way" of doing this
          do k2 = 1, nn
            M_right_der_array(k2, 1:7, 2) = I0_array(k2, 1:7, 2) * &
                                            matmul( -B_der(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 2), 6) )
            M_right_der_array(k2, 1:7, 3) = I0_array(k2, 1:7, 3) * &
                                            matmul( -B_der(1:7, k2), M_radial_monomial(lim_buffer_array(k2, 3), 6) )
          end do
          I_right_der_array(1:nn, 1:alpha_max) = matmul( M_right_der_array(1:nn, 1:7, 3), transpose(A(1:alpha_max, 1:7)) ) * &
                                                 I0_array(1:nn, 5:alpha_max + 4, 3) - &
                                                 matmul( M_right_der_array(1:nn, 1:7, 2), transpose(A(1:alpha_max, 1:7)) ) * &
                                                 I0_array(1:nn, 5:alpha_max + 4, 2)

          exp_coeff_buffer_der_array = I_left_der_array + I_right_der_array

          do j = 1, nn
            k2 = rjs_idx(j)
            exp_coeff_der(1:alpha_max, k2) = exp_coeff_der(1:alpha_max, k2) + &
                                             amplitudes(j) * exp_coeff_buffer_der_array(j, 1:alpha_max) + &
                                             amplitudes_der(j) * exp_coeff_buffer_array(j, 1:alpha_max)
          end do

          deallocate( M_left_der_array, M_right_der_array, I_left_der_array, I_right_der_array, &
                      exp_coeff_buffer_der_array, B_der )
        end if

!       deallocate( soft_weights, buffer_weights )
        deallocate( rjs, rjs_idx, atom_sigma_scaleds, s2s, atom_widths, lim_buffer_array, &
                    I0_array, g_aux_left_array, g_aux_right_array, M_left_array, M_right_array, &
                    I_left_array, I_right_array, amplitudes, amplitudes_der, B, exp_coeff_buffer_array )
      else
        deallocate( atom_widths )
      end if

!     Transform from g_alpha to g_n (the orthonormal basis) when using weights
!     We would need to allocate here again a few things!!!!!!!!!!!!
!      exp_coeff_soft_array(1:nn, 1:alpha_max) = matmul( exp_coeff_soft_array(1:nn, 1:alpha_max), W )
!      exp_coeff_buffer_array(1:nn, 1:alpha_max) = matmul( exp_coeff_buffer_array(1:nn, 1:alpha_max), W )

!      do j = 1, nn
!        k2 = rjs_idx(j)
!        exp_coeff(1:alpha_max, k2) = amplitudes(j) * soft_weights(j) * exp_coeff_soft_array(j, 1:alpha_max) + &
!                                     amplitudes(j) * buffer_weights(j) * exp_coeff_buffer_array(j, 1:alpha_max)
!      end do


!     Transform from g_alpha to g_n (the orthonormal basis)
      exp_coeff(1:alpha_max, k+1:k+n_neigh(i)) = matmul( W, exp_coeff(1:alpha_max, k+1:k+n_neigh(i)) )
      if( do_derivatives )then
        exp_coeff_der(1:alpha_max, k+1:k+n_neigh(i)) = matmul( W, exp_coeff_der(1:alpha_max, k+1:k+n_neigh(i)) )
      end if

      k = k + n_neigh(i)

!      deallocate( rjs, rjs_idx, atom_sigma_scaleds, s2s, atom_widths, lim_soft_array, soft_weights, &
!                  buffer_weights, I0_array, g_aux_left_array, g_aux_right_array, M_left_array, M_right_array, &
!                  I_left_array, I_right_array, exp_coeff_soft_array, amplitudes, amplitudes_der, B, &
!                  exp_coeff_buffer_array )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end do

!   **************** New basis ******************
!   This results from the change of variable in the
!   overlap integrals. We only need this if we want to
!   know the actual value of the expansion coefficients.
!   Since this is a global factor, once we have normalized
!   the SOAP vectors it does not have an effect anymore.
    exp_coeff = exp_coeff * dsqrt(rcut_hard_in)
    if( do_derivatives )then
      exp_coeff_der = exp_coeff_der / dsqrt(rcut_hard_in)
    end if
!   *********************************************

!   This is for debugging
    if( .false. )then
      open(10, file="radial_expansion_coefficients.dat", status="unknown")
      do i=1, size(exp_coeff, 2)
        write(10,*) exp_coeff(1:alpha_max, i)
      end do
      close(10)
      if( do_derivatives )then
        open(10, file="radial_expansion_coefficients_derivatives.dat", status="unknown")
        do i=1, size(exp_coeff, 2)
          write(10,*) exp_coeff_der(1:alpha_max, i)
        end do
        close(10)
      end if
    end if

    deallocate( A )

  return
  end subroutine
!**************************************************************************






!**************************************************************************
!
! Auxiliary function that computes a piecewise polynomial function and
! its derivatives with respect to r
!
  function g_aux(r, r0, width, piece) result(poly)
    implicit none
 
    real(dp), intent(in) :: r, r0, width
    character(*), intent(in) :: piece
    real(dp) :: x
    real(dp) :: poly(1:4)

    x = (r - r0)/width

    if ( piece == "left" )then
      poly = [ 1.d0 - 3.d0*x**2 - 2.d0*x**3, -6.d0*(x**2 + x)/width, -3.d0*(2.d0*x + 1)/width**2, -2.d0/width**3 ]
    else if ( piece == "right" )then
      poly = [ 1.d0 - 3.d0*x**2 + 2.d0*x**3, 6.d0*(x**2 - x)/width, 3.d0*(2.d0*x - 1)/width**2, 2.d0/width**3 ]
    end if

  return
  end function
!**************************************************************************

!**************************************************************************
!
! Auxiliary function that computes a piecewise polynomial function and
! its derivatives with respect to r, given an array r(:,:) of points
!
  function g_aux_array(r, r0, width, piece) result(poly)
    implicit none
 
    real(dp), intent(in) :: r(:,:), r0(:), width(:)
    character(*), intent(in) :: piece
    real(dp), dimension( 1:size(r,1), 1:size(r,2) ) :: x, x2, x3
    real(dp), dimension( 1:size(width) ) :: w, w2, w3
    real(dp), dimension( 1:size(r,1), 1:4, 1:size(r,2) ) :: poly
    integer :: i, j

    w(:) = 1.d0/width(:)

    do i = 1, size(r,2)
      x(:, i)  = (r(:, i) - r0(:)) * w
    end do
    x2 = x*x
    x3 = x2*x

    w2 = w*w
    w3 = w2*w

    if ( piece == "left" )then
      do i = 1, size(r,2)
        poly(:, 1, i) = 1.d0 - 3.d0*x2(:, i) - 2.d0*x3(:, i)
        poly(:, 2, i) = -6.d0*(x2(:, i) + x(:, i)) * w
        poly(:, 3, i) = -3.d0*(2.d0*x(:, i) + 1) * w2
        poly(:, 4, i) = -2.d0 * w3
      end do
    else if ( piece == "right" )then
      do i = 1, size(r,2)
        poly(:, 1, i) = 1.d0 - 3.d0*x2(:, i) + 2.d0*x3(:, i)
        poly(:, 2, i) = 6.d0*(x2(:, i) - x(:, i)) * w
        poly(:, 3, i) = 3.d0*(2.d0*x(:, i) - 1) * w2
        poly(:, 4, i) = 2.d0 * w3
      end do
    end if

  return
  end function
!**************************************************************************

!**************************************************************************
!
! Auxiliary function derivatives of g_aux() with respect to r0, where
!                width = width0 + width_scaling * r0
! and
!   width_scaling = 2.d0*sqrt(2.d0*log(2.d0)) * atom_sigma_scaling
!
  function g_aux_der_array(r, r0, width, width_scaling, piece) result(poly)
    implicit none
 
    real(dp), intent(in) :: r(:,:), r0(:), width(:), width_scaling
    character(*), intent(in) :: piece
    real(dp), dimension( 1:size(r,1), 1:size(r,2) ) :: x, x2, x3
    real(dp), dimension( 1:size(width) ) :: w, w2, w3, w4
    real(dp), dimension( 1:size(r,1), 1:4, 1:size(r,2) ) :: poly
    integer :: i, j

    w(:) = 1.d0/width(:)

    do i = 1, size(r,2)
      x(:, i)  = (r(:, i) - r0(:)) * w
    end do
    x2 = x*x
    x3 = x2*x

    w2 = w*w
    w3 = w2*w
    w4 = w3*w

    if ( piece == "left" )then
      do i = 1, size(r,2)
        poly(:, 1, i) = 6.d0*(x(:, i) + (width_scaling + 1.d0)*x2(:, i) + width_scaling*x3(:, i)) * w
        poly(:, 2, i) = 6.d0*(1.d0 + 2.d0*(width_scaling + 1.d0)*x(:, i) + 3.d0*width_scaling*x2(:, i)) * w2
        poly(:, 3, i) = 6.d0*(width_scaling + 1.d0 + 3.d0*width_scaling*x(:, i)) * w3
        poly(:, 4, i) = 6.d0*width_scaling * w4
      end do
    else if ( piece == "right" )then
      do i = 1, size(r,2)
        poly(:, 1, i) = 6.d0*(x(:, i) + (width_scaling - 1.d0)*x2(:, i) - width_scaling*x3(:, i)) * w
        poly(:, 2, i) = 6.d0*(1.d0 + 2.d0*(width_scaling - 1.d0)*x(:, i) - 3.d0*width_scaling*x2(:, i)) * w2
        poly(:, 3, i) = 6.d0*(width_scaling - 1.d0 - 3.d0*width_scaling*x(:, i)) * w3
        poly(:, 4, i) = -6.d0*width_scaling * w4
      end do
    end if


  return
  end function
!**************************************************************************

!**************************************************************************
!
! This subroutine returns the matrix containing the atom-independent
! coefficients for the polynomial radial basis.
! (in the notes: A, A*)
!
  subroutine get_constant_poly_coeff(alpha_max, rcut, A)
    implicit none

    integer, intent(in) :: alpha_max
    real(dp), intent(in) :: rcut
    integer :: i, j, l(1:2)
    integer, allocatable :: factor(:)
    real(dp), intent(inout) :: A(:,:)

    allocate( factor(1:alpha_max) )
    
    do i = 1, alpha_max
      A(i, 1) = rcut / N_a( rcut, i ) / dfloat(i + 3)
      factor(i) = (i + 3)
    end do

    l = shape(A)
    do j = 2, l(2)
        A(:, j) = A(:, j - 1) * rcut / dfloat(factor + j - 1) 
    end do

  end subroutine
!**************************************************************************

!**************************************************************************
!
! This subroutine returns the constant coefficients for the 
! contribution from the buffer zone, coming from the dot products
! of poly. radial basis and smoothing function (filter)
! (in the notes: B_l, B_r)
!
  subroutine get_constant_poly_filter_coeff_array(rj, width_j, r_filter, filter_width, piece, B)
    implicit none

    real(dp), intent(in) :: rj(:), width_j(:), r_filter, filter_width
    character(*), intent(in) :: piece
    integer :: i, j, k
    real(dp) :: C_filter(1:4)
    real(dp), allocatable :: col_poly(:, :, :), C_poly(:, :), rj_temp(:,:)
    real(dp), intent(out) :: B(:, :)

    allocate( col_poly(1:size(rj), 1:4, 1:1) )
    allocate( C_poly(1:7*size(rj), 1:4) )
    allocate( rj_temp(1:size(rj), 1:1) )

!   coeff. from the filter
    C_filter = g_aux(0.d0, r_filter, filter_width, "right")

!   build Toeplitz matrix a.k.a. diagonal-constant matrix
    C_poly = 0.d0

    rj_temp(:, 1) = -rj(:)

!    col_poly = g_aux_array(-reshape(rj, [size(rj), 1]), 0.d0*rj, width_j, piece)
    col_poly = g_aux_array(rj_temp, 0.d0*rj, width_j, piece)

    do i = 1, 4
      do j = 1, size(rj)
        k = (j-1)*7
        C_poly(k+i:k+i+3, i) = col_poly(j, 1:4, 1)
      end do
    end do

!    B = transpose( reshape(matmul( C_poly, C_filter ), [7, size(rj)]) )
    B = reshape(matmul( C_poly, C_filter ), [7, size(rj)])

    deallocate( col_poly, C_poly, rj_temp )

  end subroutine
!**************************************************************************

!**************************************************************************
!
! This subroutine returns the constant coefficients for the 
! contribution from the buffer zone for the derivatives, coming from
! the dot products of poly. radial basis and smoothing function (filter)
! (in the notes: partial_der(B_l), partial_der(B_r))
!
  subroutine get_constant_poly_filter_coeff_der_array(rj, width_j, width_scaling, r_filter, filter_width, piece, B)
    implicit none

    real(dp), intent(in) :: rj(:), width_j(:), r_filter, filter_width, width_scaling
    character(*), intent(in) :: piece
    integer :: i, j, k
    real(dp) :: C_filter(1:4)
    real(dp), allocatable :: col_poly(:, :, :), C_poly(:, :), rj_temp(:,:)
    real(dp), intent(out) :: B(:, :)

    allocate( col_poly(1:size(rj), 1:4, 1:1) )
    allocate( C_poly(1:7*size(rj), 1:4) )
    allocate( rj_temp(1:size(rj), 1:1) )

!   coeff. from the filter
    C_filter = g_aux(0.d0, r_filter, filter_width, "right")

!   build Toeplitz matrix a.k.a. diagonal-constant matrix
    C_poly = 0.d0

    rj_temp(:, 1) = -rj(:)

    col_poly = g_aux_der_array(rj_temp, 0.d0*rj, width_j, width_scaling, piece)

    do i = 1, 4
      do j = 1, size(rj)
        k = (j-1)*7
        C_poly(k+i:k+i+3, i) = col_poly(j, 1:4, 1)
      end do
    end do

!    B = transpose( reshape(matmul( C_poly, C_filter ), [7, size(rj)]) )
    B = reshape(matmul( C_poly, C_filter ), [7, size(rj)])

    deallocate( col_poly, C_poly, rj_temp )

  end subroutine
!**************************************************************************

!**************************************************************************
!
! This subroutine returns the radial terms coming from the
! polynomial radial basis 
! (in the notes: I0, R_hard)
!
  function M_radial_poly_array(r, a_max, rcut) result(radial_terms)
    implicit none

    integer, intent(in) :: a_max
    real(dp), intent(in) :: rcut, r(:,:)
    integer :: i, j
!   The 1st dimension of r is neighbor index and the third dimension 1:3 for the 3 limits
    real(dp), dimension(1:size(r,1), 1:a_max, 1:size(r,2)) :: radial_terms
    real(dp), dimension(1:size(r,1), 1:size(r,2)) :: r_slice

    r_slice = 1.d0 - r/rcut

    do j = 1, size(r, 2)
      radial_terms(:, 1, j) = 1.d0
      do i = 2, a_max
        radial_terms(:, i, j) = radial_terms(:, i-1, j) * r_slice(:, j)
      end do
    end do

  end function
!**************************************************************************

!**************************************************************************
!
! This function returns a vector of monomial terms (1, x, x**2, ...)
! for a given polynomial degree
! (in the notes: R, R_ext) 
!
  function radial_monomial(r, degree) result(radial_terms)
    implicit none

    integer, intent(in) :: degree
    real(dp), intent(in) :: r
    integer :: p
    real(dp), dimension(1:degree + 1) :: radial_terms 

    radial_terms = 1.d0
    do p = 2, degree + 1
      radial_terms(p) = r * radial_terms(p - 1) 
    end do

  end function
!**************************************************************************

!**************************************************************************
!
! This function returns a vector of monomial terms (1, x, x**2, ...)
! for a given polynomial degree
! (in the notes: R, R_ext) 
!
  function radial_monomial_array(r, degree) result(radial_terms)
    implicit none

    integer, intent(in) :: degree
    real(dp), intent(in) :: r(:)
    integer :: p
    real(dp), dimension(1:size(r), 1:degree + 1) :: radial_terms 

    radial_terms = 1.d0
    do p = 2, degree + 1
      radial_terms(:, p) = r(:) * radial_terms(:, p - 1) 
    end do

  end function
!**************************************************************************

!**************************************************************************
!
! This function returns the matrix of monomial terms and their
! derivatives, for a given polynomial degree:
!        ( 1  x  x**2    ...  x**degree              )
!        ( 0  1  2 * x   ...  degree * x**(degree-1) )
!        ( .   .                                     )
!        ( .       .                                 )
!        ( .             .    degree!                )
!
! (in the notes: R*) 
!
! Hardcoded for degree = 6
! Note that there is a size mismatch between radial_terms() and
! M_radial() for any other degree
!
  function M_radial_monomial(r, degree) result(M_radial)
    implicit none

    integer, intent(in) :: degree
    real(dp), intent(in) :: r
    real(dp), dimension(1:degree + 1) :: radial_terms
    real(dp), dimension(1:degree + 1, 1:degree + 1) :: M_radial

    radial_terms = radial_monomial(r, degree)

    M_radial = 0.d0
    M_radial(1:degree + 1, 1) = radial_terms(1:7)
    M_radial(2:degree + 1, 2) = radial_terms(1:6) * [1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0]
    M_radial(3:degree + 1, 3) = radial_terms(1:5) * [2.d0, 6.d0, 12.d0, 20.d0, 30.d0]
    M_radial(4:degree + 1, 4) = radial_terms(1:4) * [6.d0, 24.d0, 60.d0, 120.d0]
    M_radial(5:degree + 1, 5) = radial_terms(1:3) * [24.d0, 120.d0, 360.d0]
    M_radial(6:degree + 1, 6) = radial_terms(1:2) * [120.d0, 720.d0]
    M_radial(7:degree + 1, 7) = 720.d0

  end function
!**************************************************************************

!**************************************************************************
!
! This function returns a collection of matrices, each of them containing
! the monomial terms and their derivatives, for an array r(:) of points
! and a given polynomial degree:
!        ( 1  x  x**2    ...  x**degree              )
!        ( 0  1  2 * x   ...  degree * x**(degree-1) )
!        ( .   .                                     )
!        ( .       .                                 )
!        ( .             .    degree!                )
!
! (in the notes: R*) 
!
! Hardcoded for degree = 6
!
  function M_radial_monomial_array(r) result(M_radial_t)
    implicit none

    integer :: i
    real(dp), intent(in) :: r(:)
    real(dp), dimension(1:size(r), 1:7) :: radial_terms
    real(dp), dimension(1:size(r), 1:7, 1:7) :: M_radial
    real(dp), dimension(1:7, 1:7, 1:size(r)) :: M_radial_t
    real(dp), dimension(6) :: v1 = [1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 6.d0]
    real(dp), dimension(5) :: v2 = [2.d0, 6.d0, 12.d0, 20.d0, 30.d0]
    real(dp), dimension(4) :: v3 = [6.d0, 24.d0, 60.d0, 120.d0]
    real(dp), dimension(3) :: v4 = [24.d0, 120.d0, 360.d0]
    real(dp), dimension(2) :: v5 = [120.d0, 720.d0]

    radial_terms = radial_monomial_array(r, 6)

    M_radial = 0.d0
    M_radial(:, 1:7, 1) = radial_terms(:, 1:7)
    do i = 1, 6
      M_radial(:, i + 1, 2) = radial_terms(:, i) * v1(i)
    end do
    do i = 1, 5
      M_radial(:, i + 2, 3) = radial_terms(:, i) * v2(i)
    end do
    do i = 1, 4
      M_radial(:, i + 3, 4) = radial_terms(:, i) * v3(i)
    end do
    do i = 1, 3
      M_radial(:, i + 4, 5) = radial_terms(:, i) * v4(i)
    end do
    do i = 1, 2
      M_radial(:, i + 5, 6) = radial_terms(:, i) * v5(i)
    end do
    M_radial(:, 7, 7) = 720.d0

    do i = 1, size(r)
      M_radial_t(1:7, 1:7, i) = M_radial(i, 1:7, 1:7)
    end do

  end function
!**************************************************************************
!
!-----------------------------------------------------------------------------------
! OPERATOR FORM OF RADIAL EXP. COEFF. ENDS HERE
!-----------------------------------------------------------------------------------


end module soap_turbo_radial
