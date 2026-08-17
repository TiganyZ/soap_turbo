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

module soap_turbo_angular

! Kind parameters, declared locally on purpose: soap_turbo is a submodule with
! its own upstream and must not depend on the parent tree's kinds module.
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64



! private so a module that uses another does not inherit a second copy
  integer, parameter, private :: sp  = real32
  integer, parameter, private :: dp  = real64
  integer, parameter, private :: i32 = int32
  integer, parameter, private :: i64 = int64
  private :: real32, real64, int32, int64

! ---------------------------------------------------------------------------
! Coefficients of the Plm recursions. They are functions of (l,m) alone, and
! the original recomputed them -- including one division per element -- for
! every neighbour pair. Built once per lmax and reused.
! ---------------------------------------------------------------------------
  real(dp), allocatable, save, private :: plm_a(:), plm_b(:)      ! Plm recursion on l
  real(dp), allocatable, save, private :: der_a(:), div_a(:)      ! the two derivative recursions
  integer, save, private :: plm_lmax = -1, der_lmax = -1

  contains

!**************************************************************************
! Fills the (l,m)-only coefficients of the Plm recursion in get_plm_array,
! for every lmax up to the one asked for.
  subroutine build_plm_tables(lmax)
    implicit none
    integer, intent(in) :: lmax
    integer :: l, m, k, kmx
    if( lmax <= plm_lmax ) return
    kmx = 1 + lmax*(lmax+1)/2 + lmax
    if( allocated(plm_a) ) deallocate(plm_a, plm_b)
    allocate( plm_a(1:kmx), plm_b(1:kmx) )
    plm_a = 0.d0
    plm_b = 0.d0
    do l = 3, lmax
      do m = 0, l-2
        k = 1 + l*(l+1)/2 + m
        plm_a(k) = dfloat(2*l-1) / dfloat(l-m)
        plm_b(k) = dfloat(l-1+m) / dfloat(l-m)
      end do
    end do
    plm_lmax = lmax
  return
  end subroutine
!**************************************************************************

!**************************************************************************
! Fills the (l,m)-only coefficients of the two derivative recursions in
! get_plm_array_der.
  subroutine build_plm_der_tables(lmax)
    implicit none
    integer, intent(in) :: lmax
    integer :: l, m, k, kmx
    if( lmax <= der_lmax ) return
    kmx = 1 + lmax*(lmax+1)/2 + lmax
    if( allocated(der_a) ) deallocate(der_a, div_a)
    allocate( der_a(1:kmx), div_a(1:kmx) )
    der_a = 0.d0
    div_a = 0.d0
    do l = 0, lmax
      do m = 1, l
        k = 1 + l*(l+1)/2 + m
        der_a(k) = 0.5d0 * dfloat(l+m) * dfloat(l-m+1)
        div_a(k) = 0.5d0 * dfloat(l-m+1) * dfloat(l-m+2)
      end do
    end do
    der_lmax = lmax
  return
  end subroutine
!**************************************************************************

!**************************************************************************
! Returns an array with the Associated Legendre polynomials Plm(x), where
! Pll(x) is the highest order in the series, together with, all the lower
! order ones. l is > 0, m can take values from 0 to l and x is within the
! [-1, 1] domain plm_array is a 1D array with modified index lm -> k, where
! k = 1 + l*(l+1)/2 + m
  subroutine get_plm_array(plm_array, lmax, x)

    implicit none
    integer :: l, m, lmax, k
    real(dp) :: plm_array(:), x, sq, x2
    if(lmax < 0 .or. dabs(x) > 1.d0)then
        write(*,*) "Bad arguments for associated Legendre polynomial"
    end if
    call build_plm_tables(lmax)
    x2 = x*x
    sq = dsqrt(1.d0 - x2)
!   We need these 6 polynomials to initialize the recursion series
!   P_00(x), k = 1
    plm_array(1) = 1.d0
    if( lmax > 0 )then
!     P_10(x), k = 2
      plm_array(2) = x
!     P_11(x), k = 3
      plm_array(3) = - sq
      if( lmax > 1 )then
!       P_20(x), k = 4
        plm_array(4) = 1.5d0*x2 - 0.5d0
!       P_21(x), k = 5
        plm_array(5) = -3.d0 * x * sq
!       P_22(x), k = 6
        plm_array(6) = 3.d0 - 3.d0*x2
      else
        return
      end if
    else
      return
    end if
    if( lmax == 2 )then
      return
    else
      do l = 3, lmax
!       First we need to obtain Pl0, Pl1, up to m=l-2 with the recursion formula on l:
!       Plm(x) = ( (2l-1)*x*Pl-1m(x) - (l-1+m)*Pl-2m(x) ) / (l-m)
        do m = 0, l-2
          k = 1 + l*(l+1)/2 + m
          plm_array(k) = plm_a(k)*x*plm_array(k-l) - plm_b(k)*plm_array(k-2*l+1)
        end do
!       Now we get Pll-1 and Pll with the recursion formulas
!       Pll-1(x) = x*(2l-1)*Pl-1l-1(x)
!       Pll(x) = -(2l-1)*sqrt(1-x^2)*Pl-1l-1(x)
        k = k + 1
        plm_array(k) = x * dfloat(2*l-1) * plm_array(k-l)
        k = k + 1
        plm_array(k) = - dfloat(2*l-1) * sq * plm_array(k-l-1)
      end do
    end if
  return
  end subroutine
!**************************************************************************






!**************************************************************************
! This subroutine computes arrays of associated Legendre polynomials of the
! first kind which have been suitably modified to avoid the singularity at
! the poles (theta = 0, pi). For the derivative wrt the azimuthal angle phi
! we return -m * Plm(cos(t)) / sin(t). For m = 0, the coefficients are
! zero.
! For the derivative wrt the azimuthal angle theta, we return:
! sin(t) * d(Plm(cos(t)))/d(cos(t)).
! These tricks remove all the singularities in the gradient in spherical
! coordinates except for that at r = 0.
!
! The subroutine takes as input the Plm, which have been computed up to
! lmax+1, since that is needed in one of the recursions 
  subroutine get_plm_array_der(plm_array, lmax, x, plm_array_div_sin, plm_array_der_mul_sin)

    implicit none
    integer :: l, m, lmax, k, k_lp1_mp1, k_lp1_mm1, k_l_mp1, k_l_mm1, k_temp
    real(dp), intent(in) :: plm_array(:), x
    real(dp), intent(out) :: plm_array_div_sin(:), plm_array_der_mul_sin(:)
    real(dp) :: part1, part2

    if(lmax < 0 .or. dabs(x) > 1.d0)then
        write(*,*) "Bad arguments for associated Legendre polynomial"
    end if
    call build_plm_der_tables(lmax)

!   This is for
!   sin(t) * d(Plm(cos(t)))/d(cos(t)) = 1/2 * [ (l+m)*(l-m+1)*P_l^{m-1}(cost(t))
!                                               - P_l^{m+1}(cost(t)) ]
!   every k in 1..kmax is assigned below, so no pre-zeroing is needed here
    do l = 0, lmax
      do m = 0, l
        k = 1 + l*(l+1)/2 + m
        k_l_mp1 = k + 1
        k_l_mm1 = k - 1
!       If m = 0 then we are asking for P_l^{-1}, which is not defined. We need
!       to rewrite in terms of P_l^1:
        if( m == 0 )then
!         P_0^1 = 0
          if( l == 0 )then
            part1 = 0.d0
!         P_l^{-1} = - (l-1)!/(l+1)! * P_l^1
          else
            k_temp = 1 + l*(l+1)/2 + 1
            part1 = -0.5d0 * plm_array(k_temp)
          end if
        else
          part1 = der_a(k) * plm_array(k_l_mm1)
        end if
        if( m == l )then
          part2 = 0.d0
        else
          part2 = -0.5d0 * plm_array(k_l_mp1)
        end if
        plm_array_der_mul_sin(k) = part1 + part2
      end do
    end do

!   This is for
!   -m * Plm(cos(t))/sin(t) = 1/2 * [ (l-m+1)*(l-m+2)*P_{l+1}^{m-1}(cost(t))
!                                      - P_{l+1}^{m+1}(cost(t)) ]
!   the m = 0 entries are the only ones the loop below does not reach
    do l = 0, lmax
      plm_array_div_sin(1 + l*(l+1)/2) = 0.d0
    end do
    do l = 0, lmax
!     Note m starts at 1 here
      do m = 1, l
        k = 1 + l*(l+1)/2 + m
        k_lp1_mp1 = 1 + (l+1)*(l+2)/2 + m + 1
        k_lp1_mm1 = 1 + (l+1)*(l+2)/2 + m - 1
!       If m = 0 then we are asking for P_{l+1}^{-1}, which is not defined. We would need
!       to rewrite in terms of P_{l+1}^1, except that for m = 0 these coefficients are
!       zero anyway.
        part1 = div_a(k) * plm_array(k_lp1_mm1)
        part2 = 0.5d0 * plm_array(k_lp1_mp1)
        plm_array_div_sin(k) = part1 + part2
      end do
    end do
  return
  end subroutine
!**************************************************************************






!**************************************************************************
!
! This subroutine returns the complex conjugate of all the e^{imphi}.
! It also calls the ilexp function and multiplies its result times the
! complex exponential, returning the final prefactors which depend on l
! and m. It also lets you retrieve the product of the complex exponential
! times the radial derivative of the ilexp function.
!
! eimphi_azi_der is retained only for callers that still ask for it; the
! azimuthal derivative is i*eimphi, so get_angular_expansion_coefficients folds
! it into its own loop and does not pass this argument.
  subroutine get_eimphi_conjg(eimphi, prefl, prefm, fact_array, lmax, phi, rj, atom_sigma, &
                              scaling, do_derivatives, prefl_rad_der, eimphi_rad_der, eimphi_azi_der )
    implicit none
    real(dp), intent(in) :: rj, atom_sigma, scaling
    complex(dp) :: eimphi(:), eimphi_rad_der(:)
    complex(dp), optional :: eimphi_azi_der(:)
    real(dp) :: phi, rjbysigma, pref, cosm2, sinm2, cosm1, sinm1, cos0, sin0, cosphi2
    integer :: lmax, l, m, k
    real(dp) :: prefl(0:), fact_array(:), prefl_rad_der(0:), pref_rad_der
    complex(dp) :: prefm(0:)
    logical, intent(in) :: do_derivatives

    rjbysigma = rj/atom_sigma

!   This is fast
    call get_ilexp(prefl, fact_array, lmax, rjbysigma)
!   Complex exponential using Euler's formula and Chebyshev recursion
    cosm2 = dcos(phi)
    cosphi2 = 2.d0 * cosm2
    sinm2 = -dsin(phi)
    cosm1 = 1.d0
    sinm1 = 0.d0
    prefm(0) = (1.d0, 0.d0)
    do l = 1, lmax
      cos0 = cosphi2 * cosm1 - cosm2
      sin0 = cosphi2 * sinm1 - sinm2
      cosm2 = cosm1
      sinm2 = sinm1
      cosm1 = cos0
      sinm1 = sin0
      prefm(l) = cos0*(1.d0, 0.d0) - sin0*(0.d0, 1.d0)
    end do
    k = 1
    do l = 0, lmax
      pref = prefl(l)
      do m = 0, l
!       This is somewhat slow but not critical
        eimphi(k) = pref * prefm(m)
        k = k + 1
      end do
    end do
    if( do_derivatives )then
!     Get ilexp derivatives
      call get_ilexp_der(prefl, lmax, rj, atom_sigma, scaling, prefl_rad_der)
      k = 1
      do l = 0, lmax
        pref_rad_der = prefl_rad_der(l)
        do m = 0, l
!         This is somewhat slow but not critical
!         Radial derivative:
          eimphi_rad_der(k) = pref_rad_der * prefm(m)
          k = k + 1
        end do
      end do
!     Azimuthal angle derivative:
!     In the expression for the derivatives, there is a -i*m term multiplying the rest,
!     which we omit here, only retaining i, since we include -m as part of the Plm_array
      if( present(eimphi_azi_der) )then
        eimphi_azi_der = eimphi * dcmplx(0.d0, 1.d0)
      end if
    end if
    return
  end subroutine
!**************************************************************************






!**************************************************************************
!
! For a description see definition of the function ilexp in functions.f90
!
! The semifactorials (2i+1)!! that the recursion below needs are a function of
! lmax alone. get_angular_expansion_coefficients fills fact_array once per call
! with get_semifactorial_array; this routine only reads it.
  subroutine get_semifactorial_array(fact_array, lmax)
    implicit none
    real(dp) :: fact_array(:), fact
    integer :: i, lmax
    if( lmax > 0 )then
      fact = 1.d0
      do i = 1, lmax
        fact = fact * (2.d0*dfloat(i) + 1.d0)
        fact_array(i) = fact
      end do
    end if
  return
  end subroutine
!**************************************************************************



  subroutine get_ilexp(ilexp_array, fact_array, lmax, x)
    implicit none

    real(dp) :: ilexp_array(0:), x, x2, x4, fl, flm1, flm2, xcut = 1.d-7
    real(dp) :: inv_x2, x2pow, ratio
    integer :: l, lmax
    real(dp) :: fact_array(:)

    if( lmax < 0 )then
      write(*,*) "Bad argument (l<0) for function ilexp_double!"
      return
    end if
    x2 = x*x
    x4 = x2*x2
    if( x2 > 0.d0 )then
      inv_x2 = 1.d0/x2
    else
      inv_x2 = 0.d0
    end if
    x2pow = x2*x2                     ! x2**l, carried up with l

!   Full calculation. This is numerically unstable for small x, that's why
!   we have cases below

    if( x > 0.0d0 ) then
       flm2 = dabs( (1.d0 - dexp(-2.d0 * x2)) / 2.d0 / x2 )
       flm1 = dabs( (x2 - 1.d0 + dexp(-2.d0 * x2)*(x2+1.d0)) / 2.d0 / x4 )
    else
       flm2 = 1.0d0
       flm1 = 0.0d0
    endif

    do l = 0, lmax
      if( l == 0 )then
        if( x < xcut )then
          ilexp_array(0) = 1.d0 - x2
        else
          ilexp_array(0) = flm2
        end if
      else if( l == 1 )then
        if( x2/1000.d0 < xcut )then
          ilexp_array(1) = (x2 - x4)/ fact_array(1)
        else
          ilexp_array(1) = flm1
        end if
      else
!       x2pow is x2**l, carried up one power per step instead of being raised
!       from scratch, and the 1/x2 is formed once above
        ratio = x2pow / fact_array(l)
        if( ratio * dfloat(l) < xcut )then
          fl = ratio
        else
          fl = dabs( flm2 - (2.d0*dfloat(l) - 1.d0)*inv_x2 * flm1 )
        end if
        x2pow = x2pow * x2
        flm2 = flm1
        flm1 = fl
        ilexp_array(l) = fl
      end if
    end do
  return
  end subroutine
!**************************************************************************






!**************************************************************************
!
! This subroutine returns the radial derivatives of ilexp(rj/srj; l=0,lmax)
! If rj = 0. then the derivative is zero. The atom_sigma here is the angular
! scaled atom sigma and the scaling is the angular atom sigma scaling
!
  subroutine get_ilexp_der(ilexp_array, lmax, rj, atom_sigma, scaling, ilexp_array_der)
    implicit none

    real(dp), intent(in) :: ilexp_array(0:), rj, atom_sigma, scaling
    integer, intent(in) :: lmax
    real(dp), intent(inout) :: ilexp_array_der(0:)
    integer :: l
    real(dp) :: coeff1, coeff2

    coeff1 = 2.d0*rj/atom_sigma**2
    coeff2 = 1.d0 - scaling*rj/atom_sigma

    ilexp_array_der = 0.d0

!   Set derivative to zero if rj is (approximately) zero
    if( rj < 1.d-5 )then
      return
    end if

    ilexp_array_der(0) = coeff1 * (ilexp_array(1) - ilexp_array(0))
    do l = 1, lmax
      ilexp_array_der(l) = (-coeff1 - dfloat(2*l+2)/rj) * ilexp_array(l) + coeff1*ilexp_array(l-1)
    end do

    ilexp_array_der = ilexp_array_der * coeff2

  return
  end subroutine
!**************************************************************************
!
! This subroutine gets the angular expansion coefficients.
!
! Three things changed relative to the original scalar version, none of them
! to the mathematics:
!
!   * When derivatives are on, the Plm are computed once, to lmax+1, instead
!     of twice. The index k = 1 + l*(l+1)/2 + m does not depend on lmax, so
!     the first kmax entries of the lmax+1 array *are* the lmax array.
!
!   * The four output expressions are fused into one loop over k. They share
!     preflm*plm and the same eimphi, and the radial derivative is the value
!     again up to one factor, so the original's four whole-array expressions
!     (each with its own temporary) collapse to a handful of operations per k.
!     eimphi_azi_der is not needed at all: it is i*eimphi.
!
!   * Only the columns that are skipped get zeroed. The original zeroed all
!     four arrays up front and then overwrote nearly every column of them,
!     which for l_max = 8 is a wasted 460 MB of stores per 160k pairs.
!
  subroutine get_angular_expansion_coefficients(n_sites, n_neigh, thetas, phis, rjs, atom_sigma_in, &
                                                atom_sigma_scaling, rcut, lmax, &
                                                eimphi, preflm, plm_array, prefl, prefm, fact_array, &
                                                mask, n_species, eimphi_rad_der, do_derivatives, &
                                                prefl_rad_der, exp_coeff, exp_coeff_rad_der, &
                                                exp_coeff_azi_der, exp_coeff_pol_der )

    implicit none

    integer, intent(in) :: n_species
    integer :: lmax, kmax, n_neigh(:), n_sites, i, j, k, kk, kmax_der, i_sp
    complex(dp), intent(out) :: exp_coeff(:,:), exp_coeff_rad_der(:,:), exp_coeff_azi_der(:,:), exp_coeff_pol_der(:,:)
    real(dp) :: thetas(:), phis(:), atom_sigma_in(:), atom_sigma, atom_sigma_scaling(:), rjs(:), x, theta, phi, rj
    real(dp), intent(in) :: rcut
    real(dp) :: amplitude, rcut2, sig_fac, u, v
!   I should probably allocate and save these variables internally to minimize the number of variables that
!   need to be passed to this subroutine
    complex(dp) :: eimphi(:), eimphi_rad_der(:)
    real(dp) :: preflm(:), plm_array(:)
    real(dp) :: prefl(0:), fact_array(:), prefl_rad_der(0:)
    complex(dp) :: prefm(0:)
    logical, intent(in) :: mask(:,:), do_derivatives
    real(dp), allocatable :: plm_array_der(:), plm_array_div_sin(:), plm_array_der_mul_sin(:)
    complex(dp) :: ec, ez
    logical :: do_this

    kmax = 1 + lmax*(lmax+1)/2 + lmax
    rcut2 = rcut*rcut

!   The semifactorials the ilexp recursion needs are a function of lmax alone.
!   The original rebuilt them inside get_ilexp, once per neighbour pair.
    call get_semifactorial_array(fact_array, lmax)

    if( do_derivatives )then
      kmax_der = 1 + (lmax+1)*(lmax+2)/2 + lmax+1
      allocate( plm_array_der(1:kmax_der) )
      allocate( plm_array_div_sin(1:kmax) )
      allocate( plm_array_der_mul_sin(1:kmax) )
    end if

    k = 0
    do i = 1, n_sites
      do j = 1, n_neigh(i)
        k = k + 1
        rj = rjs(k)
        do_this = .false.
        if( rj < rcut )then
          do i_sp = 1, n_species
            if( mask(k, i_sp) )then
              do_this = .true.
              exit
            end if
          end do
        end if
        if( .not. do_this )then
!         Nothing is written for this pair, so this is the only column that
!         has to be cleared
          do kk = 1, kmax
            exp_coeff(kk, k) = (0.d0, 0.d0)
          end do
          if( do_derivatives )then
            do kk = 1, kmax
              exp_coeff_rad_der(kk, k) = (0.d0, 0.d0)
              exp_coeff_azi_der(kk, k) = (0.d0, 0.d0)
              exp_coeff_pol_der(kk, k) = (0.d0, 0.d0)
            end do
          end if
          cycle
        end if

        phi = phis(k)
        theta = thetas(k)
        x = dcos(theta)
        atom_sigma = atom_sigma_in(i_sp) + atom_sigma_scaling(i_sp)*rj
        amplitude = rcut2 / atom_sigma**2
        call get_eimphi_conjg( eimphi, prefl, prefm, fact_array, lmax, phi, rj, atom_sigma, &
                               atom_sigma_scaling(i_sp), do_derivatives, prefl_rad_der, eimphi_rad_der )

        if( do_derivatives )then
!         Plm up to lmax+1; its first kmax entries are the lmax array, which is
!         why plm_array itself is not filled again here
          call get_plm_array(plm_array_der, lmax+1, x)
          call get_plm_array_der(plm_array_der, lmax, x, plm_array_div_sin, plm_array_der_mul_sin)
          sig_fac = 2.d0 * atom_sigma_scaling(i_sp) / atom_sigma
          do kk = 1, kmax
            u  = amplitude * preflm(kk)
            v  = u * plm_array_der(kk)
            ec = v * eimphi(kk)
            exp_coeff(kk, k)         = ec
            exp_coeff_rad_der(kk, k) = v * eimphi_rad_der(kk) - sig_fac * ec
!           multiplying by i is a swap and a sign, not a complex multiply
            ez = u * plm_array_div_sin(kk) * eimphi(kk)
            exp_coeff_azi_der(kk, k) = dcmplx( -aimag(ez), dble(ez) )
            exp_coeff_pol_der(kk, k) = u * plm_array_der_mul_sin(kk) * eimphi(kk)
          end do
        else
          call get_plm_array(plm_array, lmax, x)
          do kk = 1, kmax
            exp_coeff(kk, k) = amplitude * preflm(kk) * plm_array(kk) * eimphi(kk)
          end do
        end if
      end do
    end do

    if( do_derivatives )then
      deallocate( plm_array_der, plm_array_div_sin, plm_array_der_mul_sin )
    end if

  return
  end subroutine get_angular_expansion_coefficients
!**************************************************************************




end module soap_turbo_angular
