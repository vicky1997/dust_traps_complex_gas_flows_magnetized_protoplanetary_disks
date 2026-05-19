! this module should contain all the initial conditions for dust
module initproblem

   use constants,  only: pi, third, AU
   use discstruct, only: gasmass, cs, omegaK, alpha, sigmag, xmax_gasmass
   use parameters, only: dtg, minrad0, maxrad0, r0, matdens, con2, nr, nz, ncell
   use types,      only: swarm
   use advection,  only: stokesnr
   
   implicit none
   
   private
   public   :: init_swarms, init_random_seed, mswarm, m0, nord, nmonom0, remove_swarm, add_swarm, &
                  init_virtual_box
   
   real            :: mswarm, m0, nord, nmonom0
   
   contains
   
   ! initializing the swrm array and some variables
   subroutine init_swarms(Ntot, swrm)
      implicit none
      type(swarm), dimension(:), allocatable, target   :: swrm        ! list of all swarms in the simulation
      integer, intent(in)                              :: Ntot        ! number of all swarms
      real                                             :: mdust       ! mass of the dust in the simulatated domain
      real, dimension(2)                               :: rand        ! random numbers
      real                                             :: Hg          ! pressure height scale of the gas
      integer                                          :: i
      real, parameter                                  :: s = 0.0   ! to initialize the particles density r slope (for MMSN s = -0.25)

      ! total mass of dust =  dust to gas ratio x mass of the gas               
      mdust = dtg * gasmass(minrad0*AU,maxrad0*AU,0.0)
      ! mass of one swarm
      mswarm = mdust / real(Ntot)
      ! monomer mass
      m0 = 4. * third * pi * r0**3 * matdens
      ! orders of magnitude in mass possible to cover in the simulation
      nord = (log10(mswarm/m0))
      ! initial number of monomers
      nmonom0 = mdust/m0
      
      if (.not.allocated(swrm)) allocate( swrm(Ntot) )

      ! initializing the particles
      do i = 1, Ntot
         swrm(i)%idnr = i
         swrm(i)%mass = m0 
         swrm(i)%npar = mswarm / swrm(i)%mass
         call random_number(rand)
        ! swrm(i)%rdis = minrad0 * AU + (maxrad0 - minrad0) * AU * (real(i)-0.5) / real(Ntot)   ! BE CAREFUL WITH 1D vertical column tests!!!!!!; you don't want the s then
        ! swrm(i)%rdis = minrad0 * AU + (maxrad0 - minrad0) * AU * ((s+1.)*rand(1))**(1./(s+1.)) 
         swrm(i)%rdis = (((maxrad0*AU)**(s+1.) - (minrad0*AU)**(s+1.))*rand(1) + (minrad0*AU)**(s+1.))**(1./(s+1.))
         Hg = cs(swrm(i)%rdis) / omegaK(swrm(i)%rdis)
         call random_number(rand)
         swrm(i)%zdis = Hg * sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
         swrm(i)%stnr = stokesnr(swrm(i), 0.0)
         swrm(i)%velr = 0.0
         swrm(i)%velz = 0.0
         swrm(i)%coll_f = 0
         swrm(i)%origin = 0
         swrm(I)%frag = 0
         swrm(I)%stick = 0
      enddo
          
      return
   end subroutine init_swarms

   subroutine init_virtual_box(Ntot, swrm, time, rmax_virtual, nr)
      implicit none
      type(swarm), dimension(:), allocatable, target   :: swrm        ! list of all swarms in the simulation
      type(swarm), dimension(:), allocatable, target   :: virtual_swrm
      type(swarm), dimension(:), allocatable, target   :: temp_swrm
      integer, intent(inout)                           :: Ntot        ! number of all swarms
      real, intent(inout)                              :: time        ! initial time
      real                                             :: mdust       ! mass of the dust in the simulatated domain
      real, dimension(2)                               :: rand        ! random numbers
      real                                             :: Hg, Hd      ! pressure height scale of the gas
      integer                                          :: i           ! counter
      real                                             :: x           ! random variable
      real                                             :: power, grain_size ! power law from growth barrier and grain size
      real                                             :: s = 0.          !  power index to initialize particles (s=-p+1)
      integer :: Nvirtual
      real, intent(out) :: rmax_virtual
      integer :: nr
      Nvirtual = ncell*nz
      nr = nr + 1
      mdust = Nvirtual*mswarm
      ! total mass of dust maxrad0<maxrad0+deltar
      rmax_virtual = xmax_gasmass(mdust/dtg, time, maxrad0*AU)
      open(667,file='virtualboxmass.out',position='append')
      write(*,*) "total solid mass  in the virtual box is ", mdust, "and rmax is ", rmax_virtual/AU
      close(667)
      write(*,*) "total solid mass  in the virtual box is ", mdust, "and rmax is ", rmax_virtual/AU
   
      if (.not.allocated(virtual_swrm)) allocate(virtual_swrm(Nvirtual))
   
      ! initializing the particles
      do i = 1, Nvirtual
         virtual_swrm(i)%idnr = Ntot + i
         call random_number(x)
         virtual_swrm(i)%rdis = ((rmax_virtual**(s+1.) - (maxrad0*AU)**(s+1.))*x + (maxrad0*AU)**(s+1.))**(1./(s+1.))
         virtual_swrm(i)%mass = m0   
         virtual_swrm(i)%npar = mswarm / virtual_swrm(i)%mass
         Hg = cs(swrm(i)%rdis) / omegaK(swrm(i)%rdis)
         call random_number(rand)
         virtual_swrm(i)%zdis = Hg * sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
         virtual_swrm(i)%stnr = stokesnr(swrm(i), 0.0)
         virtual_swrm(i)%velr = 0.0
         virtual_swrm(i)%velz = 0.0
         virtual_swrm(i)%coll_f = 0
         virtual_swrm(i)%origin = 1
         virtual_swrm(I)%frag = 0
         virtual_swrm(I)%stick = 0
      enddo
   
      allocate(temp_swrm(Ntot))
      temp_swrm = swrm
      deallocate(swrm)
      allocate(swrm(Ntot + Nvirtual))
   
      swrm(:Ntot) = temp_swrm
      swrm(Ntot+1:) = virtual_swrm
      deallocate(temp_swrm, virtual_swrm)
      Ntot = Ntot + Nvirtual
      return   
   end subroutine init_virtual_box

   ! add swarms with the given number
   subroutine add_swarm(swrm, Nadd, x, Ntot_glob)
      implicit none
      type(swarm), dimension(:), allocatable, target   :: swrm        ! list of all swarms in the simulation
      type(swarm), dimension(:), allocatable            :: tempswrm, tempswrm2
      integer, intent(in)                              :: Nadd        ! number of swarms to be added
      real, intent(in)                                 :: x !the location where the particles should be added
      integer, intent(inout)                           :: Ntot_glob
      real                                             :: mdust       ! mass of the dust in the simulatated domain
      real, dimension(2)                               :: rand        ! random numbers
      real                                             :: Hg, Hd      ! pressure height scale of the gas and dust, total no of particles
      integer                                          :: i, ntot
      real, parameter                                  :: s = 0.0   ! to initialize the particles density r slope (for MMSN s = -0.25)
      real                     :: stdf, stfrag, stdrift, stmax, rad, amax, rand2, q, Stnr, time
      real                     :: vfrag, gamma, dx, minrad, maxrad, rloc, dr
      real, dimension(:), allocatable  :: a_col
      real, dimension(size(swrm))      :: a_swrm
      character(len=100)               :: command
      allocate(tempswrm(Nadd))
      Ntot = size(swrm)
      !afrag = (2*sigmag(x,time)/3/pi/alpha/rho_s) 
      gamma = 11./4.
      dx = 0.005*AU
      minrad = x - dx
      maxrad = x + dx
      Hg = cs(x)/omegaK(x)

      !stfrag = (1/3./alpha(x)) * (vfrag/cs(x))**2.
      !stdf = 2 * vfrag * x * omegaK(x)/cs(x)**2/gamma
      !stdrift = (dtg/gamma) * (Hg/x)**(-2) 
      !stmax = min(stdrift, stfrag, stdf)
      !amax = 2 * sigmag(x,time) * stmax / pi / matdens
      !q = -2.5
      a_swrm = con2 * swrm(:)%mass**third
      if (maxval(a_swrm) > 5*r0) then
         !where(swrm(:)%rdis < (maxrad0 - 0.05)*AU) a_col = con2 * swrm(:)%mass**third
         a_col = pack(a_swrm,swrm(:)%rdis<(maxrad0 - 0.05)*AU)
          rloc = (maxrad0-0.05)/2
         dr = 0.05
         open(43,file='sizedistribution.dat',status='unknown',position='append')
         write(43,*) ' #mswarm '
         write(43,*) mswarm
         write(43,*) ' #rloc'
         write(43,*) rloc
         write(43,*) ' #dr'
         write(43,*) dr
         write(43,*) ' #size'
         write(43,*) a_col
         close(43)
         write(*,*)'wrting file for power law'
         write(command,*) 'python powerlaw.py'
         CALL SYSTEM(command)
         open(unit=2,file='powerlawexp.dat',action='read') !read path from file
         read(2,*) q
         close(2)

         open(unit=2,file='.dat',action='read') !read path from file
         read(2,*) amax
         close(2)

         deallocate(a_col)

         write(*,*) q, amax
         stop
      endif
      !if (stdf < stfrag) q = -2.75
      !if (stdrift < stdf) q = -2.5
      do i = 1, Nadd
         tempswrm(i)%idnr = Ntot_glob + i
         call random_number(rand)
         if (maxval(a_swrm)> 5*r0) then
            rad = (rand(1)*(amax**(4+q)-r0**(4+q))+ r0**(4+q))**(1/(4+q))
         else 
            rad = r0
         endif
         tempswrm(i)%mass = 4*third*matdens*pi*rad**3
         tempswrm(i)%npar = mswarm / tempswrm(i)%mass
         tempswrm(i)%rdis = minrad + (maxrad - minrad) * (real(i)-0.5) / real(Ntot)
         tempswrm(i)%zdis = 0.
         Stnr = stokesnr(tempswrm(i),time)
         Hd = Hg * sqrt(alpha(tempswrm(i)%rdis)/(alpha(tempswrm(i)%rdis)+Stnr))
         call random_number(rand)
         tempswrm(i)%zdis = Hd * sqrt(-2.*log(rand(1))) * cos(2.*pi*rand(2))
         tempswrm(i)%stnr = stokesnr(tempswrm(i), 0.0)
         tempswrm(i)%velr = 0.0
         tempswrm(i)%velz = 0.0
         tempswrm(i)%coll_f = 0
      enddo

      allocate(tempswrm2(Ntot))

      tempswrm2 = swrm
      deallocate(swrm)
      allocate(swrm(Ntot+Nadd))
      swrm(1:Ntot) = tempswrm2
      swrm(Ntot+1:Ntot+Nadd) = tempswrm

      deallocate(tempswrm, tempswrm2)
      Ntot_glob = Ntot_glob + Nadd
      open(63,file='feedingparticles.dat',status='unknown',position='append')
      write(63,*) Nadd, size(swrm), Ntot_glob
      close(63)    
      return
   end subroutine

   subroutine remove_swarm(swrm, Nrem)
      type(swarm), dimension(:), allocatable, target   :: swrm, tempswrm        ! list of all swarms in the simulation
      integer, intent(in)                              :: Nrem        ! number of swarms to be added
      integer                                          :: i, Nin, Nout, Nstart, Nend, Ntot, j
      real                                             :: rand

      Nin = count(swrm(:)%rdis < (minrad0)*AU)
      Ntot = size(swrm)
      allocate(tempswrm(Ntot))
      tempswrm = swrm
      deallocate(swrm)
      allocate(swrm(Ntot - Nin))
      j = 1
      do i=1,Ntot
         if(tempswrm(i)%rdis > minrad0*AU) then
            swrm(j) = tempswrm(i)
            j = j + 1
         end if
      enddo
      deallocate(tempswrm)
      open(43,file='innerparticles.dat',status='unknown',position='append')
      write(43,*) Nin, size(swrm)
      close(43)
      return
      
   end subroutine   

      ! initialize the random number generator
   subroutine init_random_seed
      implicit none
      integer                            :: i, n, clock
      integer, dimension(:), allocatable :: seed

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)

      seed = 37 * [(i - 1, i = 1, n)]
      seed = seed + clock

      call random_seed(put = seed)
      deallocate(seed)

   end subroutine init_random_seed

end
