! this module performes advection of every representative particle due to radial drift, vertical settling and turbulent
! diffusion
! in this module all the quantities should be calculated locally (at the particle location) therefore we don't need the
! grid to perform advection
module advection

   use constants,  only: smallv, pi, AU, third, AH2, mH2
   use discstruct, only: Pg, densg, omegaK, cs, diffcoefgas, vgas, ddensgdr, ddensgdz, alpha
   use types,      only: swarm
   use parameters, only: matdens, smallr, con2, vertsett, rdrift, minrad0, maxrad0

   implicit none

   private
   public      :: mc_advection, stokesnr, vel_rd, vel_vn, update_St, &
                           calculate_virtual_vr, vertical_advection

   contains

   ! routine performes radial and vertical advection due to radial drift, vertical settling and turbulent diffusion
   ! for every particle in the simulation
   subroutine mc_advection(swrm, dtime, realtime)
      implicit none
      type(swarm), dimension(:), allocatable, target :: swrm
      type(swarm), pointer                           :: particle
      real, intent(in)                               :: dtime
      real, intent(in)                               :: realtime
      real                                           :: vs, vr, velr, velv, vn
      integer                                        :: i

      !allocate(vs(size(swrm)), vr(size(swrm)), velr(size(swrm)), velv(size(swrm)), vn(size(swrm)))

      ! loop over all the particles: calculating advection velocities
      !$OMP PARALLEL DO PRIVATE(particle,vs, vr, velr, velv, vn) SCHEDULE(DYNAMIC)
      do i = 1, size(swrm)
         if (swrm(i)%rdis < smallr) cycle
         particle => swrm(i)
         particle%stnr = stokesnr(particle,realtime)
         call vel_vn(particle, vn, realtime)
         call vel_vs(particle, vs,realtime)
         call vel_ver(particle, vs, velv, dtime, realtime)
         call vel_rd(particle, vr, vn, realtime)
         call vel_rad(particle, velr, vr, dtime, realtime)

         particle%rdis = particle%rdis + velr * dtime
         if ((particle%rdis > minrad0*AU) .and. (particle%rdis<maxrad0*AU)) then
            particle%zdis = particle%zdis + velv * dtime
            particle%stnr = stokesnr(particle,realtime)
         else
            if(particle%rdis < minrad0*AU) particle%rdis = 0.5*smallr
            if(particle%rdis > maxrad0*AU) particle%rdis = 0.75*smallr
            particle%velr = 1.e-20
            particle%velz = 1.e-20
         endif
         enddo
         !$OMP END PARALLEL DO
      ! apllying the shift to particles coordinates
      !if (vertsett) then
      !   call vertical_settling(swrm, velv, dtime)
      !else
      !   call vertical_redistr(swrm, realtime)
      !endif
      !if (rdrift) then
      !   call radial_drift(swrm, velr, dtime)
      !endif

      !do i = 1, size(swrm)
      !   if (particle%rdis > smallr) particle%stnr = stokesnr(particle,realtime)
      !   if(particle%rdis < minrad0*AU) particle%rdis = 0.5*smallr
      !   if(particle%rdis > maxrad0*AU) particle%rdis = 0.75*smallr
      !enddo

      !deallocate(vs, vr, velr, velv, vn)

      return
   end subroutine mc_advection

   ! the maximum radial drift speed due to the gas pressure gradient
   subroutine vel_vn(particle, vn, realtime)
      implicit none
      type(swarm)                                    :: particle
      real, intent(in)                               :: realtime
      real                                           :: vn

      vn = 0.25 * (Pg(particle%rdis+1., particle%zdis, realtime) - Pg(particle%rdis-1.,particle%zdis, realtime)) / &
              densg(particle%rdis,particle%zdis,realtime) / omegaK(particle%rdis)

      return
   end subroutine vel_vn

   ! the complete vertical velocity
   subroutine vel_ver( particle, vs, velv, dtime, realtime)
      implicit none
      type(swarm)                                     :: particle
      !integer, intent(in)                             :: i
      real, intent(in)                                :: dtime, realtime
      real, intent(in)                                :: vs
      real                                            :: velv
      real                                            :: dz
      real                                            :: Ldiff
      real, dimension(2)                              :: rand
      real                                            :: vD1, vD2

      if (particle%rdis > smallr) then
         Ldiff = sqrt(2. * dtime * diffcoefgas(particle%rdis)  / (1. + particle%stnr**2))
         call random_number(rand)
         dz = Ldiff/(sqrt(2.*log(2.))) * sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
         vD1 = dz / dtime
         vD2 = (diffcoefgas(particle%rdis) / (1. + particle%stnr**2) / densg(particle%rdis,particle%zdis, realtime)) * &
            ddensgdz(particle%rdis,particle%zdis, realtime)
         velv = - vs + vD1 + vD2
      else
         velv = 1.e-20
      endif

      return
   end subroutine vel_ver

   ! the complete radial velocity
   subroutine vel_rad(particle, velr, vr, dtime, realtime)
      implicit none
      type(swarm)                                     :: particle
      !integer, intent(in)                             :: i
      real, intent(in)                                :: dtime, realtime
      real, intent(in)     :: vr
      real                 :: velr
      real                                            :: Ldiff
      real, dimension(2)                              :: rand
      real                                            :: dr
      real                                            :: vD1, vD2

      !if (particle%rdis > smallr) then
         Ldiff = sqrt(2. * dtime * diffcoefgas(particle%rdis) / (1. + particle%stnr**2))
         call random_number(rand)
         dr = Ldiff/(sqrt(2.*log(2.))) * sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
         vD1 = dr / dtime
         vD2 = (diffcoefgas(particle%rdis) / (1. + particle%stnr**2) / densg(particle%rdis,particle%zdis, realtime)) * &
            (ddensgdr(particle%rdis,particle%zdis, realtime))
         velr = vr + vD1 + vD2
         !particle%velr = vr
      !else
      !   velr(i) = 1.e-20
      !endif

      return
   end subroutine vel_rad

   ! performs change of the mean rdis of particles because of the radial drift
   subroutine radial_drift(swrm, velr, dtime)
      implicit none
      type(swarm), dimension(:), allocatable, target  :: swrm
      real, dimension(size(swrm))                     :: velr
      real, intent(in)                                :: dtime

      where (swrm(:)%rdis > smallr) swrm(:)%rdis = swrm(:)%rdis + velr(:) * dtime   !!!!!! vr < 0 => drift towards the star!!
      where (swrm(:)%rdis < smallr) swrm(:)%rdis = 0.99 * smallr ! artificially stops the particles right beyond the "evaporation radius"

      return
   end subroutine radial_drift

   ! performes vertical advection
   subroutine vertical_settling(swrm, velv, dtime)
      implicit none
      type(swarm), dimension(:), allocatable, target  :: swrm
      real, intent(in)                                :: dtime
      real, dimension(size(swrm))                     :: velv

      where (swrm(:)%rdis > smallr) swrm(:)%zdis = swrm(:)%zdis + velv(:) * dtime

      return
   end subroutine vertical_settling

    ! vertical redistribution of particles to a theoretical Gausian profile
   subroutine vertical_redistr(swrm, realtime)
      implicit none
      type(swarm), dimension(:), allocatable, target  :: swrm
      real, dimension(size(swrm))                     :: Hd
      real, dimension(size(swrm))                     :: rand1, rand2
      real, dimension(size(swrm))                     :: Hg
      real, intent(in)                                :: realtime
      integer                                         :: i

      do i = 1, size(swrm)
         Hg(i) = cs(swrm(i)%rdis) / omegaK(swrm(i)%rdis)
         Hd(i) = Hg(i) * (alpha(swrm(i)%rdis) / (alpha(swrm(i)%rdis)+swrm(i)%stnr))**(0.5)
      enddo
      call random_number(rand1)
      call random_number(rand2)
      where (swrm(:)%rdis > smallr) swrm(:)%zdis = Hd(:) * sqrt(-2.*log(rand1(:))) * cos(2.*pi*rand2(:))

      return
   end subroutine vertical_redistr

   ! calculates the Stokes numbers of particles locally
   real function stokesnr(particle, realtime)
      implicit none
      type(swarm)                                     :: particle
      real                                            :: lmfp
      real                                            :: rad
      real                                            :: gasdens, css
      real, intent(in)                                :: realtime

      gasdens = densg(particle%rdis,particle%zdis,realtime)
      css = cs(particle%rdis)
      lmfp = mH2 / ( gasdens * AH2 )
      rad = con2 * particle%mass**third

      if (rad > 2.25 * lmfp) then ! Stokes regime
         stokesnr = sqrt(2.*pi) * rad**2. * matdens * omegaK(particle%rdis) * AH2 / (9. * css * mH2)
      else                 ! Epstein regime
         stokesnr = rad * matdens / (sqrt(8./pi) * css * gasdens) * omegaK(particle%rdis)
      endif

      return
   end function stokesnr

   ! vertical settling velocity
   subroutine vel_vs(particle, vs, realtime)
      implicit none
      type(swarm)                                     :: particle
      real                  :: vs
      real :: realtime
      !integer, intent(in)                             :: i

      vs =particle%zdis * omegaK(particle%rdis) * particle%stnr / (1. + particle%stnr**2.)
      particle%velz = vs

      return
   end subroutine vel_vs

   ! velocity of the radial drift
   subroutine vel_rd(particle, vr, vn, realtime)
      implicit none
      type(swarm)                                     :: particle
      real, intent(in)                                :: realtime
      real                                            :: vr
      real,intent(in)     :: vn
      !integer, intent(in)                             :: i

      vr = 2. * vn / (particle%stnr + 1./particle%stnr) + &
                     vgas(particle%rdis, particle%zdis, realtime, .true.) / (1. + particle%stnr* particle%stnr)
      particle%velr = vr

      return
   end subroutine vel_rd


   subroutine update_St(swrm, realtime)
      implicit none
      type(swarm), dimension(:), allocatable, target :: swrm
      type(swarm), pointer                           :: particle
      real, intent(in)                               :: realtime
      integer                                        :: i

      !$OMP PARALLEL DO PRIVATE(particle) SCHEDULE(DYNAMIC)
      do i = 1, size(swrm)
         if (swrm(i)%rdis < smallr) cycle
         particle => swrm(i)
         particle%stnr = stokesnr(particle,realtime)
      enddo
      !$OMP END PARALLEL DO
      return
   end subroutine
   

   subroutine calculate_virtual_vr(swrm, realtime, dtime, r_VB_center, radial_vel_VB)
      implicit none
      type(swarm), dimension(:), allocatable, target :: swrm
      type(swarm), pointer                           :: particle
      real, intent(in)                               :: realtime, dtime, r_VB_center
      real, dimension(:), allocatable                :: radial_vel_VB
      real                :: vr, velr, vn
      integer                                        :: i
      real :: real_pos

      !allocate(vr(1))
      !allocate(velr(1))
      !allocate(vn(1))
      ! loop over all the particles: calculating advection velocities
      !$OMP PARALLEL DO PRIVATE(particle, vr, velr, vn, real_pos) SCHEDULE(DYNAMIC)
      do i = 1, size(swrm)

         particle => swrm(i)
         !real_pos = particle%rdis 
         !particle%rdis = r_VB_center ! all velocities calculated at the center of the virtual box
         particle%stnr = stokesnr(particle,realtime)

         call vel_vn(particle, vn ,realtime)
         call vel_rd(particle, vr, vn, realtime)
         call vel_rad(particle, velr, vr, dtime, realtime)
         radial_vel_VB(i) = velr
         !particle%rdis = real_pos ! back to particle's init pos

      enddo
      !$OMP END PARALLEL DO

      return
   end subroutine calculate_virtual_vr



   subroutine vertical_advection(swrm, realtime, dtime)
      implicit none
      type(swarm), dimension(:), allocatable, target :: swrm
      type(swarm), pointer                           :: particle
      real, intent(in)                               :: realtime, dtime
      real               :: vs, velv
      integer                                        :: i
      real :: real_pos

      !allocate(vs(1))
      !allocate(velv(1))
      ! loop over all the particles: calculating advection velocities
      !$OMP PARALLEL DO PRIVATE(particle, vs, velv) SCHEDULE(DYNAMIC)
      do i = 1, size(swrm)
         particle => swrm(i)
         particle%stnr = stokesnr(particle,realtime)
         call vel_vs(particle, vs, realtime)
         call vel_ver(particle, vs, velv, dtime, realtime)
         particle%zdis = particle%zdis + velv * dtime
      enddo
      !$OMP END PARALLEL DO

      return
   end subroutine vertical_advection

end
