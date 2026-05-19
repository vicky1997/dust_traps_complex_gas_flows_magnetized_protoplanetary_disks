!subroutines to read data from Bai2017 data for MHD wind simulation
! for more details check Bai 2017, Hu and Bai 2021
! check windvelocity.py to see how the data is adapted for our use
module Bai2017
    use constants
    use parameters
    use utils
    implicit none
    private
    public :: read_mhddata, r_bai, zh_bai, v_bai


    real, dimension(:), allocatable :: r_bai, zh_bai
    real, dimension(:,:),allocatable :: v_bai
    
    contains

    !subroutine velocity_interp(vgas, swrm)
    !end subroutine velocity_interp
    subroutine read_mhddata()
      implicit none
      real, dimension(:), allocatable :: v_cs_bai, rr_bai !values from the data file which is scale free model
      !real,dimension(:), allocatable :: z_bai, v_bai1 ! values after we scale it to our disk model
      integer :: i, stat,nr !loop variable
      real :: nrr

      open(1,file='bai2017averaged.inp',action='read')
      read(1,*) nrr
      nr = int(nrr)
      allocate(r_bai(nr),v_cs_bai(nr*nr),zh_bai(nr),rr_bai(nr*nr))
      read(1,*)
      do i=1,nr
         read(1,*) r_bai(i), zh_bai(i)
      enddo
      read(1,*)
      do i=1,nr*nr
         read(1,*) rr_bai(i), v_cs_bai(i)
      enddo
      close(1)
      !scaling the velocity and vertical height
      r_bai(:) = r_bai(:)*AU 
      !rr_bai(:) = rr_bai(:)*AU
      !do i = 1,size(r_bai)
      !   z_bai(i) = zh_bai(i)*cs(r_bai(i))/omegaK(r_bai(i))
      !enddo
      !do i =1,size(v_cs_bai)
      !   v_bai1(i) = cs(rr_bai(i))*v_cs_bai(i)
      !enddo
      allocate(v_bai(nr,nr))
      v_bai = reshape(v_cs_bai,(/nr,nr/))
      deallocate(v_cs_bai)
      !stop
      !write(*,*) r_bai(1)
      return
    end subroutine read_mhddata
   
end module Bai2017
