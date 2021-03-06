    module MassiveNu
    use constants
    use gauss_quad
    use PolyLog
    implicit none
    private

    real(dl), parameter  :: fermi_dirac_const  = 7._dl/120*const_pi**4 ! 5.68219698_dl
    !fermi_dirac_const = int q^3 F(q) dq = 7/120*pi^4
    real(dl), parameter  :: const2 = 5._dl/7._dl/const_pi**2   !0.072372274_dl

    !Steps for spline interpolation (use series outside this range)
    !terry: interpolation changed to linear spacing instead of log
    !todo check accuracy?
    integer, parameter  :: nrhopn=1400
    real(dl), parameter :: am_min = 0.3_dl
    !smallest a*m_nu to integrate distribution function rather than using series
    real(dl), parameter :: am_max = 70._dl
    !max a*m_nu to integrate

    !Actual range for using series (to avoid inaccurate ends of spline)
    real(dl), parameter :: am_minp=am_min + am_max/(nrhopn-1)*1.01_dl
    real(dl), parameter :: am_maxp=am_max*0.9_dl
    real(dl) :: xi_fac(2,4) !terry
    real(dl) :: rho_exp(3, 4) !terry
    real(dl) :: drho_exp(3, 4) !terry
    real(dl) :: P_exp(3, 4) !terry
    real(dl) :: Li357(6, 4) !terry
    real(dl) :: xi_nu(4) ! terry
    logical :: xifacOK = .false. !terry
    !integer, parameter, private :: nqmax0=80 !maximum array size of q momentum samples

    Type TNuPerturbations
        !Sample for massive neutrino momentum
        !Default settings appear to be OK for P_k accuate at 1e-3 level
        !real(dl) :: nu_q(nqmax0), nu_int_kernel(nqmax0, 4), nu_w(nqmax0) !terry
        integer nqmax !actual number of q modes evolves
        real(dl), allocatable ::  nu_q(:), nu_int_kernel(:, :), nu_w(:) !terry
    contains
    procedure :: init => TNuPerturbations_init
    end type TNuPerturbations

    Type TThermalNuBackground
        !Quantities for the neutrino background momentum distribution assuming thermal
        real(dl) dam !step in a*m
        real(dl), dimension(:,:), allocatable ::  r1,p1,dr1,dp1,ddr1 !terry
        real(dl), private :: target_rho
    contains
    procedure :: init => ThermalNuBackground_init
    procedure :: rho_P => ThermalNuBackground_rho_P
    procedure :: rho => ThermalNuBackground_rho
    procedure :: drho => ThermalNuBackground_drho
    procedure :: find_nu_mass_for_rho => ThermalNuBackground_find_nu_mass_for_rho
    end type TThermalNuBackground

    Type(TThermalNuBackground), target :: ThermalNuBackground
    class(TThermalNuBackground), pointer :: ThermalNuBack !ifort workaround

    public fermi_dirac_const,  sum_mnu_for_m1, neutrino_mass_fac, TNuPerturbations, &
        ThermalNuBackground, ThermalNuBack, nuRhoPres
    contains

    subroutine sum_mnu_for_m1(summnu,dsummnu, m1, targ, sgn)
    use constants
    real(dl), intent(in) :: m1, targ, sgn
    real(dl), intent(out) :: summnu, dsummnu
    real(dl) :: m2,m3

    m2 = sqrt(m1**2 + delta_mnu21)
    m3 = sqrt(m1**2 + sgn*delta_mnu31)
    summnu = m1 + m2 + m3 - targ
    dsummnu = m1/m2+m1/m3 + 1

    end subroutine sum_mnu_for_m1

    subroutine TNuPerturbations_init(this,Accuracy, xi_camb)
    !Set up which momenta to integrate the neutrino perturbations, depending on accuracy
    !Using three optimized momenta works very well in most cases
    class(TNuPerturbations) :: this
    real(dl), intent(in) :: Accuracy
    real(dl) :: dq,dlfdlq, q
    real(dl) :: xi_camb(4) !terry
    integer i

    this%nqmax=3
    if (Accuracy>1) this%nqmax=4
    if (Accuracy>2) this%nqmax=5
    if (Accuracy>3) this%nqmax=nint(Accuracy*10)
    !note this may well be worse than the 5 optimized points

    !We evolve evolve 4F_l/dlfdlq(i), so kernel includes dlfdlnq factor
    !Integration scheme gets (Fermi-Dirac thing)*q^n exact,for n=-4, -2..2
    !see CAMB notes
    if (allocated(this%nu_q)) deallocate(this%nu_q)
    if (allocated(this%nu_int_kernel)) deallocate(this%nu_int_kernel) ! terry
    if (allocated(this%nu_w)) deallocate(this%nu_w) ! terry
    allocate(this%nu_q(this%nqmax))
    allocate(this%nu_w(this%nqmax)) ! terry
    allocate(this%nu_int_kernel(this%nqmax, 4))

    if (this%nqmax==3) then
        !Accurate at 2e-4 level
        this%nu_q(1:3) = (/0.913201, 3.37517, 7.79184/)
        this%nu_w(1:3) = (/0.776251, 0.109268, 0.00249701/)!terry
        do i = 1, 4
          this%nu_int_kernel(1:3,i) = (/0.0687359, 3.31435, 2.29911/) !terry
        enddo
    else if (this%nqmax==4) then
        !This seems to be very accurate (limited by other numerics)
        this%nu_q(1:4) = (/0.7, 2.62814, 5.90428, 12.0/)
        this%nu_w(1:4) = (/0.747213, 0.177875, 0.0116737, 0.0000558314/)!terry
        do i = 1, 4
          this%nu_int_kernel(1:4,i) = (/0.0200251, 1.84539, 3.52736, 0.289427/)!terry
        enddo
    else if (this%nqmax==5) then
        !exact for n=-4,-2..3
        !This seems to be very accurate (limited by other numerics)
        this%nu_q(1:5) = (/0.583165, 2.0, 4.0, 7.26582, 13.0/)
        this%nu_w(1:5) = (/0.681808, 0.222159, 0.0454693, 0.00294858, 0.0000177609/)!terry
        do i = 1, 4
          this%nu_int_kernel(1:5,i) = (/0.0081201, 0.689407, 2.8063, 2.05156, 0.126817/)!terry
        enddo
    else
        stop 'not supported yet, nqmax too big, AccuracyBoost>3' !terry
        dq = (12 + this%nqmax/5)/real(this%nqmax)
        do i=1,this%nqmax
            q=(i-0.5d0)*dq
            this%nu_q(i) = q
            dlfdlq=-q/(1._dl+exp(-q))
            this%nu_int_kernel(i,1:4)=dq*q**3/(exp(q)+1._dl) * (-0.25_dl*dlfdlq) !now evolve 4F_l/dlfdlq(i)
        end do
    end if
    if (this%nqmax>=3 .and. this%nqmax<=5) then !terry
      do i = 1, this%nqmax
        this%nu_int_kernel(i, :) = this%nu_w(i) * exp(this%nu_q(i)) * this%nu_q(i)**4 * exp(this%nu_q(i)-xi_camb(:)) / &
          (1.d0 + exp(this%nu_q(i)-xi_camb(:)))**2 / 4 !terry
        this%nu_int_kernel(i, :) = this%nu_int_kernel(i,:) + this%nu_w(i) * exp(this%nu_q(i)) * this%nu_q(i)**4 *&
          exp(this%nu_q(i)+xi_camb(:)) / (1.d0 + exp(this%nu_q(i)+xi_camb(:)))**2 / 4
        this%nu_int_kernel(i,:) = this%nu_int_kernel(i,:)/2
      enddo
    endif !terry
    this%nu_int_kernel=this%nu_int_kernel/fermi_dirac_const

    end subroutine TNuPerturbations_init

    subroutine ThermalNuBackground_init(this, xi_camb)
      use constants ! terry
    class(TThermalNuBackground) :: this
    !  Initialize interpolation tables for massive neutrinos.
    !  Use cubic splines interpolation of log rhonu and pnu vs. log a*m.
    integer i, j
    real(dl) am, rhonu,pnu
    real(dl) spline_data(nrhopn), exp_int, ex
    real(dl), intent(in) :: xi_camb(4) !Kenny: edited terry

    call setPolyLogPrec(10)
    do i = 1, 4
      !terms for 0-th and second order expansion of rho, P, drho
	    xi_fac(1,i) = (1._dl+30._dl*(xi_camb(i)**2._dl)/7._dl/(const_pi**2._dl)+15._dl*(xi_camb(i)**4._dl)/7._dl/(const_pi**4._dl)) !terry
	    xi_fac(2,i) = (1._dl + 3._dl / const_pi / const_pi * xi_camb(i)**2) !terry
      !calculate terms for higher order expansion of rho, P, drho
      !tested by am=0.3, xi=1.2, errors for rho and P are about 10^-6, error for drho is about 10^-4
      !the three coefficients are for m^4*log(m), m^4 and m^5 respectively
      !alpha(a) used in the expansion is taken to be 0.001, seems accurate enough
      !remember to normalize by 7pi^4/120=fermi_dirac_const
      !todo drho terms are different from the original even for xi=0, bug?
      
      !rho terms:
      !1/16
      ![1-4*log(2a)]/64 + a*exp(x)/8/(1+exp(x))^2 - 1/16 * int(1/p*[1/(1+exp(1/p+x))+1/(1+exp(1/p-x))], 0, 1/a)
      !-2*exp(x)/15/(1+exp(x))^2
      !P terms: note that factor of 3 is divided at the end, now just multiply them by 3
      !-1/16 * 3
      !([-7+12*log(2a)]/192 - a*exp(x)/8/(1+exp(x))^2 + 1/16 * int(1/p*[1/(1+exp(1/p+x))+1/(1+exp(1/p-x))], 0, 1/a)) * 3
      !8exp(x)/45/(1+exp(x))^2 * 3
      !drho terms: 1/4
      ![1-2*log(2a)]/8 + a*exp(x)/2/(1+exp(x))^2 - 1/4 * int(1/p*[1/(1+exp(1/p+x))+1/(1+exp(1/p-x))], 0, 1/a)
      !-2*exp(x)/3/(1+exp(x))^2
      !note: that integral is the same for the 3 cases
      ex = exp(xi_camb(i))
      Li357(:, i) = get_PolyLog_357(ex)
      ex = ex / (1.d0 + ex)**2
      exp_int = get_exp_int(xi_camb(i))
      rho_exp(1, i) = 1.d0 / 16.d0
      rho_exp(2, i) = (1.d0 - 4.d0 * log(2.d0 * 1.d-3)) / 64.d0 + 1.d-3 * ex / 8 - 1.d0 / 16.d0 * exp_int
      rho_exp(3, i) = -2.d0 / 15.d0 * ex
      P_exp(1, i) = -1.d0 / 16.d0
      P_exp(2, i) = (-7.d0 + 12.d0 * log(2.d0 * 1.d-3)) / 192.d0 - 1.d-3 * ex / 8 + 1.d0 / 16.d0 * exp_int
      P_exp(3, i) = 8.d0 / 45.d0 * ex 
      drho_exp(1, i) = 1.d0 / 4.d0
      drho_exp(2, i) = (1.d0 - 2.d0 * log(2.d0 * 1.d-3)) / 8.d0 + 1.d-3 * ex / 2 - 1.d0 / 4.d0 * exp_int
      drho_exp(3, i) = -2.d0 / 3.d0 * ex
    enddo
    rho_exp = rho_exp / fermi_dirac_const
    !for P: note that factor of 3 is divided at the end, now just multiply them by 3
    P_exp = P_exp / fermi_dirac_const * 3
    drho_exp = drho_exp / fermi_dirac_const
    xi_nu = xi_camb
    xifacOK = .true.

    !terry bug fix
    !if (allocated(this%r1)) return
    ThermalNuBack => ThermalNuBackground !ifort bug workaround

    if (.not. allocated(this%r1)) then
      allocate(this%r1(nrhopn, 4),this%p1(nrhopn, 4),this%dr1(nrhopn, 4),this%dp1(nrhopn, 4),this%ddr1(nrhopn, 4))
    endif
    this%dam=(am_max-am_min)/(nrhopn-1)

    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), &
    !$OMP& PRIVATE(am,rhonu,pnu)
    do i=1,nrhopn
        am=am_min + (i-1)*this%dam
        do j = 1, 4 !terry
          call nuRhoPres(am,rhonu,pnu, xi_camb(j)) !terry
          this%r1(i,j)=rhonu !terry
          this%p1(i,j)=pnu !terry
        enddo
    end do
    !$OMP END PARALLEL DO

    do i = 1, 4 !terry
      call splini(spline_data,nrhopn) !terry
      call splder(this%r1(:,i),this%dr1(:,i),nrhopn,spline_data) !terry
      call splder(this%p1(:,i),this%dp1(:,i),nrhopn,spline_data) !terry
      call splder(this%dr1(:,i),this%ddr1(:,i),nrhopn,spline_data) !terry
    enddo !terry

    !terry debug
    !print *, "xinu", xi_nu(1:3)
    !do j = 1, 3
    !  print *, "thermal neutrino rho p drho", j
    !  do i = -10, 10
    !    if (i==0) then
    !      print *, "middle min"
    !    endif
    !    am = am_minp*(1+i*1.d-2)
    !    call ThermalNuBackground_rho_P(this,am,rhonu,pnu, j)
    !    print *, rhonu, pnu, ThermalNuBackground_drho(this, am, 1.d0, j)
    !  enddo
    !  print *, "++++++++++++++++++"
    !  do i = -10, 10
    !    if (i==0) then
    !      print *, "middle max"
    !    endif
    !    am = am_maxp*(1+i*1.d-2)
    !    call ThermalNuBackground_rho_P(this,am,rhonu,pnu, j)
    !    print *, rhonu, pnu, ThermalNuBackground_drho(this, am, 1.d0, j)
    !  enddo
    !  print *, "++++++++++++++++++"
    !  do i = -10, 10
    !    if (i==0) then
    !      print *, "middle 0.01"
    !    endif
    !    am = 0.01d0*(1+i*1.d-2)
    !    call ThermalNuBackground_rho_P(this,am,rhonu,pnu, j)
    !    print *, rhonu, pnu, ThermalNuBackground_drho(this, am, 1.d0, j)
    !  enddo
    !enddo

    end subroutine ThermalNuBackground_init

    subroutine nuRhoPres(am,rhonu,pnu, xi_camb) !terry
    !  Compute the density and pressure of one eigenstate of massive neutrinos,
    !  in units of the mean density of one flavor of massless neutrinos.
    real(dl),  parameter :: qmax=30._dl
    integer, parameter :: nq=100
    real(dl) dum1(nq+1),dum2(nq+1)
    real(dl), intent(in) :: am
    real(dl), intent(out) ::  rhonu,pnu
    integer i
    real(dl) q,aq,v,aqdn,adq
    real(dl), intent(in) :: xi_camb !terry
!   Kenny: simply take average?
    real(dl) aqdn_anti,adq_anti              !Kenny: edited
    real(dl) dum1_anti(nq+1),dum2_anti(nq+1) !Kenny: edited
    real(dl) rhonu_anti, pnu_anti            !Kenny: edited

    !  q is the comoving momentum in units of k_B*T_nu0/c.
    !  Integrate up to qmax and then use asymptotic expansion for remainder.
    adq=qmax/nq
    dum1(1)=0._dl
    dum2(1)=0._dl
    dum1_anti(1)=0._dl !Kenny: edited
    dum2_anti(1)=0._dl !Kenny: edited
    do  i=1,nq
        q=i*adq
        aq=am/q
        v=1._dl/sqrt(1._dl+aq*aq)
!        aqdn=adq*q*q*q/(exp(q)+1._dl) !Kenny: edited
        aqdn=adq*q*q*q/(exp(-xi_camb)*exp(q)+1._dl) !Kenny: edited
        aqdn_anti=adq*q*q*q/(exp(xi_camb)*exp(q)+1._dl) !Kenny: edited
        dum1(i+1)=aqdn/v
        dum2(i+1)=aqdn*v
        dum1_anti(i+1)=aqdn_anti/v !Kenny: edited
        dum2_anti(i+1)=aqdn_anti*v !Kenny: edited
    end do
    call splint(dum1,rhonu,nq+1)
    call splint(dum2,pnu,nq+1)
    call splint(dum1_anti,rhonu_anti,nq+1)   !Kenny: edited
    call splint(dum2_anti,pnu_anti,nq+1)     !Kenny: edited
    !  Apply asymptotic corrrection for q>qmax and normalize by relativistic
    !  energy density.
    rhonu=(rhonu+dum1(nq+1)/adq)/fermi_dirac_const
    pnu=(pnu+dum2(nq+1)/adq)/fermi_dirac_const/3._dl
    rhonu_anti=(rhonu_anti+dum1_anti(nq+1)/adq)/fermi_dirac_const       !Kenny: edited
    pnu_anti=(pnu_anti+dum2_anti(nq+1)/adq)/fermi_dirac_const/3._dl     !Kenny: edited

    rhonu=(rhonu+rhonu_anti)/2._dl  !Kenny: edited   !Kenny: simply take average?
    pnu=(pnu+pnu_anti)/2._dl        !Kenny: edited   !Kenny: simply take average?
    end subroutine nuRhoPres

    subroutine ThermalNuBackground_rho_P(this,am,rhonu,pnu, nu_i) !terry
    class(TThermalNuBackground) :: this
    real(dl), intent(in) :: am
    integer, intent(in) :: nu_i !terry
    real(dl), intent(out) :: rhonu, pnu
    real(dl) d, logam, am2
    real(dl) :: approx_p1, approx_p2, x2, xi
    integer i
    !  Compute massive neutrino density and pressure in units of the mean
    !  density of one eigenstate of massless neutrinos.  Use cubic splines to
    !  interpolate from a table. Accuracy generally better than 1e-5.

    xi = xi_nu(nu_i)
    if (am <= am_minp) then
        if (am< 0.01_dl) then
            rhonu=xi_fac(1,nu_i) + const2*am**2 * xi_fac(2,nu_i)
            pnu=(2 * xi_fac(1,nu_i) - rhonu)/3._dl
        else
            !Higher order expansion result less obvious, Appendix A of arXiv:0911.2714
            am2=am**2
            logam = log(am)
            !rhonu = 1+am2*(const2+am2*(.1099926669d-1*logam-.3492416767d-2-.5866275571d-2*am))
            !pnu = (1+am2*(-const2+am2*(-.3299780009d-1*logam-.5219952794d-3+.2346510229d-1*am)))/3
            rhonu = xi_fac(1,nu_i)+am2*(const2*xi_fac(2,nu_i)+am2*(rho_exp(1, nu_i)*logam+rho_exp(2, nu_i)+rho_exp(3, nu_i)*am))
            pnu = (xi_fac(1,nu_i)+am2*(-const2*xi_fac(2,nu_i)+am2*(P_exp(1, nu_i)*logam+P_exp(2, nu_i)+P_exp(3, nu_i)*am)))/3
        end if
        return
    !else if (am >= am_maxp .and. xi/=0.d0) then
    else if (am >= am_maxp) then
      !x2=xi**2
      !approx_p1 = zeta3 + 2.d0/3.d0*log(2.d0)*x2 + x2**2/72.d0 - x2**3/4320.d0 + x2**4/120960.d0
      !approx_p2 = 15.d0*zeta5 + 6.d0*zeta3*x2 + 2.d0/3.d0*log(2.d0)*x2**2 + x2**3/180.d0 - x2**4/20160.d0
      !rhonu = 3/(2*fermi_dirac_const)*(approx_p1*am + (approx_p2)/2/am)
      !approx_p1 = zeta5 + 2.d0/5.d0*zeta3*x2 + 2.d0/45.d0*log(2.d0)*x2**2 + x2**3/2700.d0 - x2**4/302400.d0
      !approx_p2 = 63.d0/4*zeta7 + 7.5d0*x2*zeta5 + 0.5d0*zeta3*x2**2 + 1.d0/45.d0*log(2.d0)*x2**3 + x2**4/10080.d0
      !pnu = 900._dl/120._dl/fermi_dirac_const*(approx_p1-approx_p2/am**2)/am
      rhonu = 1.d0/fermi_dirac_const*(-sum(Li357(1:2,nu_i))*am + (-6.d0*sum(Li357(3:4,nu_i)) + 45._dl*sum(Li357(5:6,nu_i))/am**2)/am)
      pnu = 1._dl/fermi_dirac_const*(-4.d0*sum(Li357(3:4,nu_i))+60._dl*sum(Li357(5:6,nu_i))/am**2)/am
      return
    else if (am >= am_maxp) then
        !Simple series solution (expanded in 1/(a*m))
        rhonu = 3/(2*fermi_dirac_const)*(zeta3*am + ((15*zeta5)/2 - 945._dl/16*zeta7/am**2)/am)
        pnu = 900._dl/120._dl/fermi_dirac_const*(zeta5-63._dl/4*Zeta7/am**2)/am
        return
    end if

    d=(am-am_min)/this%dam+1._dl
    i=int(d)
    d=d-i

    !  Cubic spline interpolation.
    rhonu=this%r1(i, nu_i)+d*(this%dr1(i, nu_i)+d*(3._dl*(this%r1(i+1, nu_i)-this%r1(i, nu_i))-2._dl*this%dr1(i, nu_i) &
        -this%dr1(i+1, nu_i)+d*(this%dr1(i, nu_i)+this%dr1(i+1, nu_i)+2._dl*(this%r1(i, nu_i)-this%r1(i+1, nu_i)))))
    pnu=this%p1(i, nu_i)+d*(this%dp1(i, nu_i)+d*(3._dl*(this%p1(i+1, nu_i)-this%p1(i, nu_i))-2._dl*this%dp1(i, nu_i) &
        -this%dp1(i+1, nu_i)+d*(this%dp1(i, nu_i)+this%dp1(i+1, nu_i)+2._dl*(this%p1(i, nu_i)-this%p1(i+1, nu_i)))))

    end subroutine ThermalNuBackground_rho_P

    subroutine ThermalNuBackground_rho(this,am,rhonu, nu_i)
    class(TThermalNuBackground) :: this
    real(dl), intent(in) :: am
    real(dl), intent(out) :: rhonu
    real(dl) d, am2
    integer i, nu_i
    real(dl) :: approx_p1, approx_p2, x2, xi

    !  Compute massive neutrino density in units of the mean
    !  density of one eigenstate of massless neutrinos.  Use series solutions or
    !  cubic splines to interpolate from a table.

    xi = xi_nu(nu_i)
    if (am <= am_minp) then
        if (am < 0.01_dl) then
            rhonu = xi_fac(1,nu_i) + const2*am**2 * xi_fac(2,nu_i) !terry
        else
            am2=am**2
            rhonu = xi_fac(1,nu_i)+am2*(const2*xi_fac(2,nu_i)+am2*(rho_exp(1, nu_i)*log(am)+rho_exp(2, nu_i)+rho_exp(3, nu_i)*am))
        end if
        return
    !else if (am >= am_maxp .and. xi/=0.d0) then
    else if (am >= am_maxp) then
      !x2=xi**2
      !approx_p1 = zeta3 + 2.d0/3.d0*log(2.d0)*x2 + x2**2/72.d0 - x2**3/4320.d0 + x2**4/120960.d0
      !approx_p2 = 15.d0*zeta5 + 6.d0*zeta3*x2 + 2.d0/3.d0*log(2.d0)*x2**2 + x2**3/180.d0 - x2**4/20160.d0
      !rhonu = 3/(2*fermi_dirac_const)*(approx_p1*am + (approx_p2)/2/am)
      rhonu = 1.d0/fermi_dirac_const*(-sum(Li357(1:2,nu_i))*am + (-6.d0*sum(Li357(3:4,nu_i)) + 45._dl*sum(Li357(5:6,nu_i))/am**2)/am)
      return
    else if (am >= am_maxp) then
        rhonu = 3/(2*fermi_dirac_const)*(zeta3*am + ((15*zeta5)/2 - 945._dl/16*zeta7/am**2)/am)
        return
    end if

    d=(am-am_min)/this%dam+1._dl
    i=int(d)
    d=d-i

    !  Cubic spline interpolation.
    rhonu=this%r1(i, nu_i)+d*(this%dr1(i, nu_i)+d*(3._dl*(this%r1(i+1, nu_i)-this%r1(i, nu_i))-2._dl*this%dr1(i, nu_i) &
        -this%dr1(i+1, nu_i)+d*(this%dr1(i, nu_i)+this%dr1(i+1, nu_i)+2._dl*(this%r1(i, nu_i)-this%r1(i+1, nu_i)))))

    end subroutine ThermalNuBackground_rho

    function rho_err(this, nu_mass, nu_i)
    class(TThermalNuBackground) :: this
    real(dl) rho_err, nu_mass, rhonu
    integer :: nu_i

    call this%rho(nu_mass, rhonu, nu_i)
    rho_err = rhonu - this%target_rho

    end function rho_err

    function ThermalNuBackground_find_nu_mass_for_rho(this,rho, nu_i) result(nu_mass)
    !  Get eigenstate mass given input density (rho is neutrino density in units of one massless)
    !  nu_mass=m_n*c**2/(k_B*T_nu0).
    !  Get number density n of neutrinos from
    !  rho_massless/n = int q^3/(1+e^q) / int q^2/(1+e^q)=7/180 pi^4/Zeta(3)
    !  then m = Omega_nu/N_nu rho_crit /n if non-relativistic
    use MathUtils
    use config
    class(TThermalNuBackground) :: this
    real(dl), intent(in) :: rho
    real(dl) nu_mass, rhonu, rhonu1, delta
    real(dl) fzero
    real(dl) massless_limit
    integer iflag, nu_i

    massless_limit = xi_fac(1,nu_i) !terry
    if (rho <= 1.001_dl * massless_limit) then
        !energy density all accounted for by massless result
        nu_mass=0
    else
        !Get mass assuming fully non-relativistic
        nu_mass=fermi_dirac_const/(-sum(Li357(1:2,nu_i)))*rho !terry

        if (nu_mass>4) then
            !  perturbative correction for velocity when nearly non-relativistic
            !  Error due to velocity < 1e-5 for mnu~0.06 but can easily correct (assuming non-relativistic today)
            !  Note that python does not propagate mnu to omnuh2 consistently to the same accuracy; but this makes
            !  fortran more internally consistent between input and computed Omega_nu h^2

            !Make perturbative correction for the tiny error due to the neutrino velocity
            call this%rho(nu_mass, rhonu, nu_i)
            call this%rho(nu_mass*0.9, rhonu1, nu_i)
            delta = rhonu - rho
            nu_mass = nu_mass*(1 + delta/((rhonu1 - rhonu)/0.1) )
        else
            !Directly solve to avoid issues with perturbative result when no longer very relativistic
            this%target_rho = rho
            !call brentq(this,rho_err,0._dl,nu_mass,0.01_dl,nu_mass,fzero,iflag)
            call brentq(this,rho_e,0._dl,nu_mass,0.01_dl,nu_mass,fzero,iflag)
            if (iflag/=0) call GlobalError('find_nu_mass_for_rho failed to find neutrino mass')
        end if
    end if

    contains
      function rho_e(this, nu_mass)
        class(TThermalNuBackground) :: this
        real(dl) rho_e, nu_mass

        rho_e = rho_err(this, nu_mass, nu_i)
      end function rho_e
    end function ThermalNuBackground_find_nu_mass_for_rho

    function ThermalNuBackground_drho(this,am,adotoa, nu_i) result (rhonudot)
    !  Compute the time derivative of the mean density in massive neutrinos
    class(TThermalNuBackground) :: this
    real(dl) adotoa,rhonudot
    real(dl) d, am2
    real(dl), intent(IN) :: am
    integer i, nu_i
    real(dl) :: approx_p1, approx_p2, x2, xi

    xi = xi_nu(nu_i)
    if (am< am_minp) then
        !rhonudot = 2*const2*am**2*adotoa
        !rhonudot = rhonudot * xi_fac(2,nu_i) !terry	
        am2 = am**2
        !rhonudot=  am**2*(2*const2 +am2*(0.1759882671_dl*log(am)+.3211546524d-1 -0.1466568893_dl*am))*adotoa
        rhonudot=  am**2*(2*const2*xi_fac(2,nu_i) +am2*(drho_exp(1,nu_i)*log(am)+drho_exp(2,nu_i)+drho_exp(3,nu_i)*am))*adotoa
        rhonudot=  am**2*(2*const2 +am2*(0.0439971_dl*log(am) - 0.00297043_dl -0.0293314_dl*am))*adotoa
    !else if (am>am_maxp .and. xi/=0.d0) then
    else if (am>am_maxp) then
      !x2=xi**2
      !approx_p1 = zeta3 + 2.d0/3.d0*log(2.d0)*x2 + x2**2/72.d0 - x2**3/4320.d0 + x2**4/120960.d0
      !approx_p2 = 15.d0*zeta5 + 6.d0*zeta3*x2 + 2.d0/3.d0*log(2.d0)*x2**2 + x2**3/180.d0 - x2**4/20160.d0
      rhonudot = 1.d0/fermi_dirac_const*(-sum(Li357(1:2,nu_i))*am + ( 6.d0*sum(Li357(3:4,nu_i)) - 135.d0*sum(Li357(5:6,nu_i))/am**2)/am)*adotoa
    else if (am>am_maxp) then
        rhonudot = 3/(2*fermi_dirac_const)*(zeta3*am +( -(15*zeta5)/2 + 2835._dl/16*zeta7/am**2)/am)*adotoa
    else
        d=(am-am_min)/this%dam+1._dl
        i=int(d)
        d=d-i
        !  Cubic spline interpolation for rhonudot.
        rhonudot=this%dr1(i, nu_i)+d*(this%ddr1(i, nu_i)+d*(3._dl*(this%dr1(i+1, nu_i)-this%dr1(i, nu_i)) &
            -2._dl*this%ddr1(i, nu_i)-this%ddr1(i+1, nu_i)+d*(this%ddr1(i, nu_i)+this%ddr1(i+1, nu_i) &
            +2._dl*(this%dr1(i, nu_i)-this%dr1(i+1, nu_i)))))

        rhonudot=adotoa*rhonudot/this%dam*am
    end if

    end function ThermalNuBackground_drho

    function get_exp_int(x)
      real(dl) :: x, get_exp_int
      intent(in) :: x
      !print *, 'test func', f(0.000889347923469268663244d0)
      get_exp_int = gauss_leg_int_divide(f, 0.d0, 1.d3, 100)
    contains
      function f(p)
        real(dl) :: f, p, pinv
        intent(in) :: p
        pinv = 1.d0 / p
        !prevents overflow, when p is small enough f is 0
        !value taken for 10^300
        if (pinv > 690.d0) then
          f = 0.d0
          return
        endif
        f = pinv * (1.d0 / (1.d0 + exp(pinv + x)) + 1.d0 / (1.d0 + exp(pinv - x)))
      end function
    end function

    function get_PolyLog_357(ex)
      real(dl) :: get_PolyLog_357(6), ex, exinv
      intent(in) :: ex
      exinv = 1.d0 / ex
      get_PolyLog_357 = [polyLogLi(3, -ex), polyLogLi(3, -exinv), polyLogLi(5, -ex), &
        polyLogLi(5, -exinv), polyLogLi(7, -ex), polyLogLi(7, -exinv)]
    end function get_PolyLog_357

    end module MassiveNu
