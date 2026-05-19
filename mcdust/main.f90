     ! This code performes 2D simulation of dust evolution in protoplanetary disk.
! The gas disk is treated 100% analytically (discstruct module)
! The dust is treated as representative particles (RPs) undergoing advection (advection module)
! as well as collisions performed with Monte Carlo algorithm (collisions module)
! To perform collisions the RPs are binned using an adaptive grid (grid module)
!
! citation: Drążkowska, Windmark & Dullemond (2013) A&A 556, A37
!
! Author: Joanna Drążkowska, Heidelberg University, ZAH, ITA
! Albert-Ueberle-Str. 2, 69120 Heidelberg, Germany
! mailto: asiadrazkowska@gmail.com
!
program main
   use constants
   use advection,    only: mc_advection, calculate_virtual_vr, vertical_advection, update_St
   use collisions,   only: mc_collisions
   use grid,         only: make_grid, deallocate_grid
   use initproblem,  only: init_swarms, init_random_seed, m0, mswarm, nord, remove_swarm, add_swarm, init_virtual_box
   use discstruct,   only: alpha, cs, omegaK, sigmag, densg, Pg, vgas, gasmass, dlogPg, Mdot
   use parameters,   only: read_parameters, Ntot, nz, nr, dtime, fout, tend, smallr, restart, restime, minrad0, &
                           maxrad0, matdens, r0, dtg, db_data, mhdwinds, collswitch, con2, dumpfile
   use output,       only: write_output, read_restart
   use timestep,     only: time_step
   use types
   use hdf5
   use hdf5output, only: hdf5_file_write, hdf5_file_t, hdf5_file_read, read_ncolls, write_ncolls, write_dumpfile
   use Bai2017, only: read_mhddata

   implicit none

   ! the array of the representative particles (swarms) is declared here:
   type(swarm), dimension(:), allocatable, target                  :: swrm        ! list of all swarms in the simulation
   type(grid_container)               :: g
   type(list_of_swarms), dimension(:,:), allocatable, target       :: bin         ! swarms binned into cells
   type(list_of_swarms), dimension(:), allocatable                 :: rbin        ! swarms binned into radial zones
   type(hdf5_file_t)                                               :: file
   real                       :: total
   real(kind=4), dimension(2) :: elapsed
   ! initialize similar quantities for virtual box 
   type(grid_container)               :: g_VB
   type(list_of_swarms), dimension(:,:), allocatable, target       :: bin_VB
   type(list_of_swarms), dimension(:), allocatable                 :: rbin_VB        ! swarms binned into radial zones
   type(swarm), dimension(:), allocatable, target :: real_swrm ! swarms inside rmax, swarms in
   type(swarm), dimension(:), allocatable, target :: virtual_swrm ! swarms in virtual box
   type(swarm), dimension(:), allocatable, target :: temp_swrm ! temporary list for reallocations
   integer :: Nreal, Nvirtual
   real :: r_VB_center, deltar_VB ! center of the virtual box to compute radial velocities, and width of virtual box
   real :: r_VB_max ! outer edge of virtual box
   real, dimension(:), allocatable :: radial_vel_VB ! radial velocities of particles in the virtual box
   integer, dimension(:), allocatable :: n_arr ! number of copies of virtual particles per time step
   real :: n_bar ! helper to compute n_arr
   real :: ti ! helper to compute time required to leave virtual box
   integer :: k ! particle counter
   !real :: totmass_VB
   real :: rmax_virtual

   integer             :: i, j, iter
   real                :: time = 0.0          ! physical time
   real                :: timeofnextout = 0.0 ! time of next output
   real                :: resdt = 0.1*year               ! resulting physical time step
   integer             :: nout = 0            ! number of the next output
   real                :: totmass             ! total mass of dust beyond evaporation line
   character(len=100)  :: ctrl_file           ! parameter file
   real                :: mdust               ! mass of dust ! TODO: can be merged with totmass?
   integer, dimension(:,:), allocatable :: ncolls, ncolls_VB
   real                :: mdust_add
   integer             :: Nadd, Ntot_glob
   real                :: rand
   integer             :: lastcolidx
   !real, dimension(:), allocatable  :: mgrid  ! mass grid 
 
   ! random number generator initialization
   call init_random_seed

   ! reading parameters and initializing the simulation
   call get_command_argument(1, ctrl_file)
   write(*,*) 'MCDust v0.1'
   write(*,*) '------------------------------------------------------------------'
   write(*,*) 'Reading parameters...'
   call read_parameters(ctrl_file)
   write(*,*) '------------------------------------------------------------------'
   write(*,*) 'Initializing representative bodies...'
   if (restart) then
      write(*,*) ' Reading restart...'
      !call read_restart(Ntot, swrm)
      call hdf5_file_read(swrm, nout, mswarm, time, Ntot, resdt)
      do i = 1,size(swrm)
         swrm(i)%npar = mswarm/swrm(i)%mass
      enddo
      write(*,*) time/year, nout, mswarm, Ntot, resdt/year
      write(*,*) '  restart read!'
      !time = restime * year
      !nout = nint(time/dtime)
      timeofnextout = time + dtime
      !mdust = dtg * gasmass(minrad0*AU,maxrad0*AU,0.0)  ! TODO
      !mswarm = mdust / real(Ntot)
      mdust = mswarm * real(size(swrm))
      m0 = 4. * third * pi * r0**3 * matdens
      nord = (log10(mswarm/m0))
      Nvirtual = count(swrm(:)%rdis>maxrad0*AU)
      Nreal = size(swrm) - Nvirtual
      allocate(real_swrm(Nreal), virtual_swrm(Nvirtual))
      allocate(n_arr(Nvirtual))
      write(*,*)Nreal, Nvirtual
      ! they are radially ordered as they are saved before advection
      real_swrm(:) = swrm(:Nreal)
      virtual_swrm(:) = swrm(Nreal+1:)

      allocate(radial_vel_VB(Nvirtual))
      r_VB_max = maxval(virtual_swrm(:)%rdis)

      deltar_VB =  r_VB_max-(maxrad0*AU)
      r_VB_center = 0.5*(r_VB_max + (maxrad0*AU))
      write(*,*) 'values:', maxval(real_swrm(:)%rdis)/AU, r_VB_center/AU
      
   else
      call init_swarms(Ntot,swrm)

      call init_virtual_box(Ntot, swrm, time, rmax_virtual, nr)

      Nvirtual = count(swrm(:)%rdis>(maxrad0*AU))
      Nreal = Ntot - Nvirtual
      allocate(real_swrm(Nreal), virtual_swrm(Nvirtual))
      allocate(n_arr(Nvirtual))
      ! they are radially ordered as they are saved before advection
      real_swrm(:) = swrm(:Nreal)
      virtual_swrm(:) = swrm(Nreal+1:)
   
      allocate(radial_vel_VB(Nvirtual))
      r_VB_max = maxval(virtual_swrm(:)%rdis)
   
      deltar_VB =  r_VB_max-(maxrad0*AU)
      r_VB_center = 0.5*(r_VB_max + (maxrad0*AU))
   endif
   write(*,*) 'succeed'

   if(mhdwinds) then
      write(*,*)'Reading gas velocity data from Bai 2017 for mhd winds'
      call read_mhddata()
   endif
   write(*,*) 'Initial disk mass: ', gasmass(0.1*AU,maxrad0*AU,0.0)/Msun

   write(*,*) ' Making grid for the first time... starting with main grid'
   !write(*,*) maxval(real_swrm(:)%rdis)/AU, maxval(virtual_swrm(:)%rdis)/AU
   !write(*,*) minval(real_swrm(:)%rdis)/AU, minval(virtual_swrm(:)%rdis)/AU
   call make_grid(real_swrm, bin, rbin, nr, nz, smallr, totmass, ncolls, g)
   write(*,*) ' Making virtual grid'
   call make_grid(virtual_swrm, bin_VB, rbin_VB, 1, nz, smallr, totmass, ncolls_VB,g_VB)
   write(*,*) '  grid done'
   
   !allocate(ncolls(nr,nz))
   !if (restart) then
   !   call read_ncolls(ncolls, ncolls_VB)
   !else
   ncolls(:,:) = 1
   ncolls_VB (:,:) = 1
   !endif
   Ntot_glob = size(swrm)

   write(*,*) 'going into the main loop...'

   iter = 0
   ! ------------------- MAIN LOOP -------------------------------------------------------------------------------------
   do while (time < tend)

      
       ! determining the time step
      if (iter == 0) then 
         if (restart .eqv. .false.) resdt = 1./omegaK(minval(swrm(:)%rdis))
      else
         call time_step(bin, ncolls, timeofnextout-time, resdt, g)      
      !else
      !   call time_step(bin, ncolls, timeofnextout-time, resdt, g)
      end if
      !write(*,*) resdt/year
      write(*,*) ' Performing advection: timestep',resdt/year,'yrs'
      if (db_data) then
         open(23,file='timestep.dat',status='unknown',position='append')
         write(23,*) time/year,resdt/year
         close(23)
      endif
      !stop
      if(resdt==0.)  then
         write(*,*) 'timestep 0.'
         stop
      endif
      ! producing output
      if (modulo(iter,fout) == 0 .or. time>=timeofnextout) then
         deallocate(swrm)
         allocate(swrm(Nreal + Nvirtual))
         swrm(:Nreal) = real_swrm
         swrm(Nreal + 1:) = virtual_swrm
         call update_St(swrm, time)
         !call write_output(swrm, nout)
         call hdf5_file_write(file, swrm, time, 'create', nout, mswarm, Ntot, resdt)
         write(*,*) 'Time: ', time/year, 'produced output: ',nout
         open(23,file='timesout.dat',status='unknown',position='append')
         write(23,*) 'time: ', time/year, 'produced output: ',nout
         close(23)
         call write_ncolls(ncolls, ncolls_VB)
         timeofnextout = time+dtime
         nout = nout + 1
      endif

      iter = iter + 1

     

      ! writing max mass value for each timestep for bug fixes
      if(db_data) then
         open(123,file='mmax.dat',position='append')
         write(123,*) time/year, maxval(swrm(:)%mass)
         close(123)
      endif



     !write(*,*) 'near advection woth timestep', resdt/year, ncolls
      !if(iter==3) stop
      ! performing advection
      call mc_advection(real_swrm, resdt, time)
      write(*,*) '  advection done'
   

      call calculate_virtual_vr(virtual_swrm, time, resdt, r_VB_center, radial_vel_VB)
      open(669,file='Vbvr.out',position='append')
      write(669,*) radial_vel_VB
      close(669)
      ! vertical displacement occurs same way as for real particles
      call vertical_advection(virtual_swrm, time, resdt)
      do i=1, Nvirtual
         if (radial_vel_VB(i)<0.) then
            ti = (virtual_swrm(i)%rdis - (maxrad0*AU))/abs(radial_vel_VB(i))
            if (ti<resdt) then
               n_bar = (resdt-ti)*abs(radial_vel_VB(i))/deltar_VB + 1.
               n_arr(i) = int(n_bar) 
               virtual_swrm(i)%rdis = r_VB_max - (n_bar - int(n_bar))*deltar_VB
            else
               n_arr(i) = 0.
               virtual_swrm(i)%rdis = virtual_swrm(i)%rdis + radial_vel_VB(i)*resdt 
            endif
         else
            n_arr(i) = 0. ! it cannot create particles on this way
            ti = (r_VB_max-virtual_swrm(i)%rdis)/radial_vel_VB(i)
            if (ti<resdt) then
               n_bar = (resdt-ti)*radial_vel_VB(i)/deltar_VB + 1.
               virtual_swrm(i)%rdis = (maxrad0*AU) + (n_bar - int(n_bar))*deltar_VB
            else
               
               virtual_swrm(i)%rdis = virtual_swrm(i)%rdis + radial_vel_VB(i)*resdt 
            endif
         endif
      enddo


      if (sum(n_arr)>0) then 
         allocate(temp_swrm(Nreal))
         temp_swrm = real_swrm
         deallocate(real_swrm)
         allocate(real_swrm(Nreal + sum(n_arr)))
         real_swrm(:Nreal) = temp_swrm 
         deallocate(temp_swrm)
         k = 0
         do i=1, Nvirtual
            do j=1, n_arr(i)
               k = k+1
               real_swrm(Nreal+k) = virtual_swrm(i)
               real_swrm(Nreal+k)%idnr = Ntot + k
               real_swrm(Nreal+k)%rdis = maxrad0*AU - (j-1)*deltar_VB - (r_VB_max - virtual_swrm(i)%rdis)
            enddo
         enddo
         Ntot = Ntot + k
         Nreal = Nreal + k
         open(667,file='sumnarray.out',position='append')
         write(667,*) 'test: ', k , ' and ', sum(n_arr) , ' should be equal'
         close(667)
         open(668,file='narr.out',position='append')
         write(668,*) n_arr
         close(668)
         write(*,*) 'test: ', k , ' and ', sum(n_arr) , ' should be equal'
         open(669,file='flux.out',position='append')
         write(669,*) sum(n_arr)*mswarm, resdt/year
         close(669)
      endif
      
      ! testing
      do i = 1, Nvirtual
         if ((virtual_swrm(i)%rdis < (maxrad0*AU)) .or. (virtual_swrm(i)%rdis>r_VB_max)) then
            write(*,*) 'Problem with particle outside ', virtual_swrm(i)%rdis/AU, r_VB_max/AU
            write(*,*) ' its radial vel is ', radial_vel_VB(i)
            stop
         endif

      enddo


      ! removing old grid and building new one
      call deallocate_grid(g)
      if (allocated(bin))  deallocate(bin)
      if (allocated(rbin)) deallocate(rbin)
      if (allocated(ncolls)) deallocate(ncolls)
      call deallocate_grid(g_VB)
      if (allocated(bin_VB))  deallocate(bin_VB)
      if (allocated(rbin_VB)) deallocate(rbin_VB)
      if (allocated(ncolls_VB)) deallocate(ncolls_VB)

      write(*,*) '    Making grid...'
      call make_grid(real_swrm, bin, rbin, nr, nz, smallr,totmass, ncolls,g)
      write(*,*) ' making virtual grid'
      call make_grid(virtual_swrm, bin_VB, rbin_VB, 1, nz, smallr,totmass, ncolls_VB, g_VB)

      write(*,*) '     grid done'
      if(collswitch) then
         ! performing collisions
         write(*,*) '   Performing collisions...'

         !$OMP PARALLEL DO PRIVATE(i,j) SCHEDULE(DYNAMIC)
         do i = 1, size(g%rce)
            do j = 1, size(g%zce,dim=2)
               if (.not.allocated(bin(i, j)%p)) cycle
               !write(*,*) '    entering zone',i,j,'including ',size(bin(i, j)%p),' rbs','.....'
               call mc_collisions(i, j, bin, real_swrm, resdt, time, ncolls(i,j),g)
               !write(*,*)'out of collisons subroutine'
            enddo
         enddo
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO PRIVATE(j) SCHEDULE(DYNAMIC)
         do j = 1, size(g_VB%zce,dim=2)

            if (.not.allocated(bin_VB(1, j)%p)) cycle
            !write(*,*) '    entering virtual zone',j,'including ',size(bin_VB(1, j)%p),' rbs','.....'
            call mc_collisions(1, j, bin_VB, virtual_swrm, resdt, time, &
                               ncolls_VB(1,j), g_VB)

         enddo
         !$OMP END PARALLEL DO

         write(*,*) '    collisions done!'
      end if
      if (dumpfile .and. modulo(iter,100)==0) then
         deallocate(swrm)
         allocate(swrm(Nreal + Nvirtual))
         swrm(:Nreal) = real_swrm
         swrm(Nreal + 1:) = virtual_swrm
         call update_St(swrm, time)
         !call write_output(swrm, nout)
         call write_dumpfile(file, swrm, time, nout, mswarm, Ntot, resdt)
      endif 
      time = time + resdt
      ! checking max sizes of real and virtual swarm and writing to file

      !lastcolidx = minloc(real_swrm(:)rdis,mask=real_swrm(:) .GT. 22.5)

      open(666,file='maxgrainsize.dat',position='append')
      write(666,*) time/year, maxval(con2 * real_swrm(:)%mass**third ), maxval(con2 * virtual_swrm(:)%mass**third)
      close(666)

   enddo
   ! ---- END OF THE MAIN LOOP -----------------------------------------------------------------------------------------
   !nout = nout+1
   write(*,*) 'time: ', time/year, 'produced output: ',nout
   open(23,file='timesout.dat',status='unknown',position='append')
   !call write_output(swrm, nout)
   write(23,*) 'time: ', time/year, 'produced output: ',nout
   close(23)
   deallocate(swrm)
   allocate(swrm(Nreal + Nvirtual))
   swrm(:Nreal) = real_swrm
   swrm(Nreal + 1:) = virtual_swrm
   call update_St(swrm, time)
   call hdf5_file_write(file, swrm, time, 'create', nout, mswarm, Ntot, resdt)
   call write_ncolls(ncolls, ncolls_VB)
   
   
   deallocate(bin)
   deallocate(rbin)
   deallocate(bin_VB)
   deallocate(swrm)
   deallocate(rbin_VB, ncolls_VB)
   deallocate(radial_vel_VB)
   deallocate(virtual_swrm, real_swrm, n_arr)

   write(*,*) '------------------------------------------------------------------'
   write(*,*) 'tend exceeded, finishing simulation...'

   ! this causes problems when used with intel compilers, so just remove it in case you want to use one
   total = etime(elapsed)
   write(*,*) 'Elapsed time [s]: ', total, ' user:', elapsed(1), ' system:', elapsed(2)

end
