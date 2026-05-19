! put your gas disk model here
! x is distance from the star (in cm)
module discstruct
    use constants, only: AU, kB, mH2, Msun, year
    use parameters, only: alphaMRI, sigmag0, temperature, eta, mhdwinds
    use utils
    use Bai2017, only: r_bai, zh_bai, v_bai
    implicit none

    private
    public   :: alpha, cs, omegaK, sigmag, densg, Pg, vgas, ddensgdz, ddensgdr, gasmass, dlogPg, diffcoefgas, z_exp, Mdot &
                                ,xmax_gasmass


    !real, parameter             :: alphaMRI = 1.e-3                     ! Shakura-Sunyaev's turbulence parameter; MRI turbulence ONLY!
    real, parameter             :: q       = 0.5                        ! T propto r^-q
    real, parameter             :: p       = 1.                         ! gas surface density propto r^-p
    !real, parameter             :: T0      = 280.0                      ! temperature at 1 AU [K]
    !real, parameter             :: cs0     = sqrt(kB * temperature / mH2)        ! sound speed at 1 AU
    !real, parameter             :: sigmag1AU = 800.                     ! gas surface density at 1AU [g cm^-2]
    real, parameter             :: dM      = 1.e-9 * Msun / year        ! gas accretion rate
    contains

    ! Shakura-Sunyaev's turbulence parameter
    real function alpha(x)
        implicit none
        real, intent(in)  :: x

        alpha = alphaMRI ! this just means constant alpha

        return
    end function

    ! speed of the sound in gas
    real function cs(x)
        implicit none
        real, intent(in)  :: x
        real :: cs0

        cs0     = sqrt(kB * temperature / mH2)
        cs = cs0 * (x/(1.*AU))**(-0.5 * q)

        return
    end function

    ! Keplerian frequency
    real function omegaK(x)
        use constants, only: Ggrav, Msun
        implicit none
        real, intent(in)  :: x

        omegaK = sqrt(Ggrav * Msun / x**3)

        return
    end function

    ! gas surface density
    real function sigmag(x, time)
        use constants, only: smallv
        implicit none
        real, intent(in) :: x
        real, intent(in) :: time ! in seconds

        sigmag = sigmag0 * (x/AU)**(-1.*p)
        sigmag = max(sigmag, smallv)
        
        return
    end function
    !calculate the exponential in density
    real function z_exp(x,z, time)
        implicit none
        real, intent(in) :: x,z, time
        real :: Hg
        real, parameter :: min_exp = 3.7266e-6
         Hg = cs(x) / omegaK(x) ! gas disk scaleheight
        z_exp = exp(-0.5 * (z / Hg)**2)
        z_exp = max(z_exp, min_exp)
        return
    end function

    
    ! gas volume density
    real function densg(x,z, time)
        use constants, only: pi
        implicit none
        real, intent(in)  :: x, z
        real  :: Hg
        real, intent(in)  :: time

        Hg = cs(x) / omegaK(x) ! gas disk scaleheight
        densg = (sigmag(x, time) / (sqrt(2.*pi) * Hg)) * z_exp(x,z,time)
        !densg = max(densg, 1.e-20)
        return
    end function

    ! gas pressure
    real function Pg(x, z,time)
        implicit none
        real, intent(in)  :: x, z, time

        Pg = densg(x,z, time) * cs(x)**2.

        return
    end function

    ! radial gas velocity
    real function vgas(x, z, time, vgaswind)
        use constants, only: AU
        use Bai2017
        implicit none
        real, intent(in)  :: x, z, time
        logical, intent(in) :: vgaswind
        real :: zh, vgasalpha

        

        if(vgaswind) then
            if(mhdwinds) then
                           
                zh = z/(cs(x)/omegaK(x))
                vgasalpha = -3. * diffcoefgas(x) / 2. / x / cs(x)
                if (abs(zh) .le. 3) then

                !if(abs(zh) .le. maxval(zh_bai) .and. x/AU .ge. minval(r_bai) .and. x/AU .le. maxval(r_bai)) then
                !write(*,*) r_bai(1), zh_bai(1)
                    vgas = interp2d(r_bai/AU,zh_bai,v_bai,x/AU,zh,.false.,vgasalpha)*cs(x)
                !write(*,*) 'vgas:', vgas, 'r(AU):', x/AU, 'z (Hg)', zh, 'z(AU)',z/AU
                !stop
                else 
                    vgas = vgasalpha * cs(x)
                endif 
            else 
                vgas = -3. * diffcoefgas(x) / 2. / x
            endif 
        else
            !write(*,*) 'vgas 0'
            !vgas = -1. * dM / (2. * pi * sigmag(x, time) * x)
            vgas = 0.0
        endif
        return
    end function

    ! turbulent diffusion coefficient
    real function diffcoefgas(x)
        implicit none
        real, intent(in)  :: x

        diffcoefgas = alpha(x) * cs(x)**2 / omegaK(x)

        return
    end function

    ! partial derivative of gas density with respect to vertical height
    real function ddensgdz(x,z, time)
        use constants, only: pi
        implicit none
        real, intent(in)  :: x, z, time
        real              :: Hg

        Hg = cs(x) / omegaK(x)
        !ddensgdz = -z * sigmag(x, time) * exp(-0.5 * (z / Hg)**2) / (sqrt(2.*pi) * Hg**3)
        ddensgdz = -z * densg(x,z,time)/ Hg**2
        return
    end function

    ! partial derivative of gas density with respect to r
    real function ddensgdr(x,z, time)
        implicit none
        real, intent(in)  :: x, z, time
        real              :: dr

        dr = 0.000001 * x
        ddensgdr = ( densg(x+dr,z, time) - densg(x-dr,z, time) ) / (2. * dr)

        return
    end function

    real function dlogPg(x, z, time)
        implicit none
        real, intent(in)  :: x, time,z
        real :: Deltar = 1.e4

        dlogPg = 0.5 * x * (Pg(x+Deltar, z, time) - Pg(x-Deltar,z,time)) / Pg(x,z,time) / Deltar

        return
    end function

    ! mass of the gas included between xmin and xmax
    real function gasmass(xmin,xmax,time)
        use constants, only: pi
        use parameters, only: dtg
        implicit none
        real, intent(in)  :: xmin,xmax, time
        real              :: x
        integer           :: i
        integer, parameter:: N = 10000
        real              :: dx

        dx = (xmax - xmin) / real(N)

        gasmass = 0.0
        x = xmin + 0.5 * dx
        do i = 1, N-1
            gasmass = gasmass + sigmag(x,time) * 2.*pi*x*dx
            x = x + dx
        enddo

        return
    end function
    ! calculate feeding flux at a given radial point
    real function Mdot(x, time)
        use constants, only: pi
        use parameters, only: dtg, vfrag

        implicit none
        real, intent(in) :: x, time
        real :: velr, stfrag, stdrift, stdf, stokesnr, gamma, Hg, vn
        gamma = 11./4.
        Hg = cs(x)/omegaK(x)
        stfrag = (1/3./alpha(x)) * (vfrag/cs(x))**2.
        stdf = 2 * vfrag * x * omegaK(x)/cs(x)**2/gamma
        stdrift = (dtg/gamma) * (Hg/x)**(-2) 
        stokesnr = min(stdrift, stfrag, stdf)
        vn = 0.25 * (Pg(x+1., 0., time) - Pg(x-1.,0., time)) / &
              densg(x,0.,time) / omegaK(x)
        velr = 2. * vn * stokesnr / (stokesnr**2. + 1.) + &
        vgas(x, 0., time, .false.) / (1. + stokesnr * stokesnr)

        Mdot = - 2 * pi * x * velr * dtg * sigmag(x,time)
        return
end function

real function xmax_gasmass(value, time, xmin)
      use constants, only: pi
      implicit none
      real, intent(in)  :: xmin, value, time
      real              :: dx = 0.001*AU ! I don't think higher resolution is necessary
      real :: gasmass
      gasmass = 0.
      xmax_gasmass = xmin + 0.5*dx
      do while (.True.)
         gasmass = gasmass + sigmag(xmax_gasmass,time) * 2.*pi*xmax_gasmass*dx
         xmax_gasmass = xmax_gasmass + dx
         if (gasmass > value) exit
      enddo

      return
   end function


end module discstruct
